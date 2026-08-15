from pathlib import Path
import csv
import json
import networkx as nx

gfa_path = Path(snakemake.input.gfa)
copies_path = Path(snakemake.input.copies)
out_gfa = Path(snakemake.output.gfa)
out_copies = Path(snakemake.output.copies)
out_summary = Path(snakemake.output.summary)
out_gfa.parent.mkdir(parents=True, exist_ok=True)

mode = str(snakemake.params.mode)
radius = int(snakemake.params.radius)
min_nodes = int(snakemake.params.min_nodes)
max_nodes = int(snakemake.params.max_nodes)
repeat_threshold = float(snakemake.params.repeat_copy_threshold)

copy_rows = {}
copy_order = []
with copies_path.open(newline="") as handle:
    for r in csv.DictReader(handle, delimiter="\t"):
        copy_rows[r["node"]] = r
        copy_order.append(r["node"])

segments = {}
segment_order = []
links = []
other = []
G = nx.Graph()

with gfa_path.open() as handle:
    for raw in handle:
        line = raw.rstrip("\n")
        f = line.split("\t")
        if not f:
            continue
        if f[0] == "S":
            segments[f[1]] = line
            segment_order.append(f[1])
            G.add_node(f[1])
        elif f[0] == "L" and len(f) >= 5:
            links.append((f[1], f[3], line))
            G.add_edge(f[1], f[3])
        elif f[0] in {"H"}:
            other.append(line)

if not segments:
    raise ValueError("No GFA segments found.")

if mode == "whole" or len(segments) <= max_nodes:
    selected = set(segment_order)
    reason = "whole_graph"
else:
    repeats = {
        n for n in segment_order
        if float(copy_rows.get(n, {"copy_number": 0})["copy_number"]) >= repeat_threshold
    }
    branch_nodes = {n for n in G.nodes if G.degree(n) > 2}
    self_loops = {n for n in G.nodes if G.has_edge(n, n)}
    seeds = branch_nodes | repeats | self_loops

    candidate_sets = []
    for seed in seeds:
        lengths = nx.single_source_shortest_path_length(G, seed, cutoff=radius)
        candidate_sets.append(set(lengths))

    # Merge overlapping neighborhoods so a tangle is considered as one region.
    merged = []
    for s in candidate_sets:
        hit = [i for i, m in enumerate(merged) if m & s]
        if not hit:
            merged.append(set(s))
        else:
            union = set(s)
            for i in reversed(hit):
                union |= merged.pop(i)
            merged.append(union)

    def score(nodes):
        sub = G.subgraph(nodes)
        branch = sum(1 for n in sub.nodes if sub.degree(n) > 2)
        repeat = sum(
            1 for n in sub.nodes
            if float(copy_rows.get(n, {"copy_number": 0})["copy_number"])
            >= repeat_threshold
        )
        cycle_rank = max(0, sub.number_of_edges() - sub.number_of_nodes()
                         + nx.number_connected_components(sub))
        weight_sum = sum(
            float(copy_rows.get(n, {"copy_number": 0})["copy_number"])
            for n in sub.nodes
        )
        return 10 * branch + 5 * cycle_rank + 2 * repeat + weight_sum

    valid = [s for s in merged if len(s) >= min_nodes]
    if valid:
        selected = max(valid, key=score)
        reason = "branch_repeat_neighborhood"
    else:
        # Fallback: center a connected window on the most graph-complex node.
        center = max(
            G.nodes,
            key=lambda n: (
                G.degree(n),
                float(copy_rows.get(n, {"copy_number": 0})["copy_number"])
            )
        )
        selected = set(nx.single_source_shortest_path_length(
            G, center, cutoff=max(radius, 1)
        ))
        reason = "fallback_neighborhood"

    # Connected trimming to keep the HUBO intentionally small.
    if len(selected) > max_nodes:
        sub = G.subgraph(selected)
        center = max(
            sub.nodes,
            key=lambda n: (
                sub.degree(n),
                float(copy_rows.get(n, {"copy_number": 0})["copy_number"])
            )
        )
        keep = {center}
        frontier = [center]
        while frontier and len(keep) < max_nodes:
            u = frontier.pop(0)
            neigh = sorted(
                [v for v in sub.neighbors(u) if v not in keep],
                key=lambda v: (
                    sub.degree(v),
                    float(copy_rows.get(v, {"copy_number": 0})["copy_number"])
                ),
                reverse=True,
            )
            for v in neigh:
                if len(keep) >= max_nodes:
                    break
                keep.add(v)
                frontier.append(v)
        selected = keep
        reason += "_trimmed"

# Preserve original segment ordering.
selected_order = [n for n in segment_order if n in selected]

with out_gfa.open("w") as out:
    for line in other:
        out.write(line + "\n")
    for n in selected_order:
        out.write(segments[n] + "\n")
    for u, v, line in links:
        if u in selected and v in selected:
            out.write(line + "\n")

with out_copies.open("w") as out:
    out.write("node\tlength\tdepth\tcopy_number\n")
    for n in selected_order:
        r = copy_rows[n]
        out.write(
            f'{n}\t{r["length"]}\t{r["depth"]}\t{r["copy_number"]}\n'
        )

sub = G.subgraph(selected)
summary = {
    "mode": mode,
    "selection_reason": reason,
    "nodes": selected_order,
    "node_count": len(selected_order),
    "edge_count": sub.number_of_edges(),
    "branch_nodes": [n for n in selected_order if sub.degree(n) > 2],
    "cycle_rank": max(
        0,
        sub.number_of_edges() - sub.number_of_nodes()
        + nx.number_connected_components(sub)
    ) if sub.number_of_nodes() else 0,
}
with out_summary.open("w") as out:
    json.dump(summary, out, indent=2)
