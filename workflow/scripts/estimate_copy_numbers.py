from pathlib import Path
import csv
import json
import math
import networkx as nx

gfa = Path(snakemake.input.gfa)
depths_tsv = Path(snakemake.input.depths)
out_tsv = Path(snakemake.output.copies)
out_json = Path(snakemake.output.baseline)
out_tsv.parent.mkdir(parents=True, exist_ok=True)

mode = str(snakemake.params.mode)
min_depth_initial = float(snakemake.params.min_depth)
depth_div = float(snakemake.params.depth_div)
offset = float(snakemake.params.offset)
min_copy = float(snakemake.params.min_copy)
max_copy = float(snakemake.params.max_copy)
fit_step = float(snakemake.params.fit_step)

rows = []
with depths_tsv.open(newline="") as handle:
    for r in csv.DictReader(handle, delimiter="\t"):
        rows.append({
            "node": r["node"],
            "length": int(r["length"]),
            "covered_bases": float(r["covered_bases"]),
            "depth": float(r["depth"]),
        })

if not rows:
    raise ValueError("No node depths were produced.")

G = nx.Graph()
for r in rows:
    G.add_node(r["node"])

with gfa.open() as handle:
    for raw in handle:
        f = raw.rstrip("\n").split("\t")
        if not f:
            continue
        if f[0] == "L" and len(f) >= 5:
            G.add_edge(f[1], f[3])

cycle_nodes = set()
try:
    for cycle in nx.cycle_basis(G):
        cycle_nodes.update(cycle)
except Exception:
    pass
self_loop_nodes = {u for u, v in nx.selfloop_edges(G) if u == v}

def weighted_mean(items):
    den = sum(r["length"] for r in items)
    if den <= 0:
        return 0.0
    return sum(r["depth"] * r["length"] for r in items) / den

# qpg's estimator deliberately avoids loops/self-loops when estimating
# single-copy depth, because repeats can inflate their coverage.
eligible = [
    r for r in rows
    if r["node"] not in cycle_nodes
    and r["node"] not in self_loop_nodes
    and r["depth"] > min_depth_initial
]
if not eligible:
    eligible = [r for r in rows if r["depth"] > min_depth_initial]
if not eligible:
    raise ValueError(
        "All nodes have near-zero depth; mapping/copy-number estimation cannot continue."
    )

if mode == "median":
    vals = sorted(r["depth"] for r in eligible)
    mid = len(vals) // 2
    baseline = (
        vals[mid] if len(vals) % 2
        else (vals[mid-1] + vals[mid]) / 2.0
    )
    fit_delta = None
else:
    # Iterative weighted depth estimate, following the logic in
    # qpg/tag_gfa_copy_numbers.pl.
    threshold = min_depth_initial
    baseline = weighted_mean([r for r in eligible if r["depth"] > threshold])
    for _ in range(10):
        if baseline <= 0:
            break
        threshold = baseline / depth_div
        used = [r for r in eligible if r["depth"] > threshold]
        if not used:
            break
        new_baseline = weighted_mean(used)
        if abs(new_baseline - baseline) < 1e-12:
            baseline = new_baseline
            break
        baseline = new_baseline

    # Fit the fundamental one-copy depth by looking for a depth D for which
    # observed depths are near integer multiples of D. This mirrors the
    # search in qpg/tag_gfa_copy_numbers.pl.
    lo = baseline / 1.3
    hi = baseline * 1.5
    if fit_step <= 0:
        fit_step = max(baseline / 100.0, 1e-6)

    best = baseline
    best_delta = math.inf
    candidate = lo
    while candidate <= hi + fit_step * 0.5:
        if candidate > 0:
            delta = 0.0
            for r in rows:
                d = r["depth"]
                if d > baseline / 4.0:
                    multiple = max(1, int(d / candidate + 0.5))
                    diff = abs(d - multiple * candidate) / candidate
                    delta += diff * r["length"]
            if delta < best_delta:
                best_delta = delta
                best = candidate
        candidate += fit_step
    baseline = best
    fit_delta = best_delta

if baseline <= 0:
    raise ValueError(f"Invalid inferred one-copy depth: {baseline}")

for r in rows:
    copy_number = r["depth"] / baseline + offset
    copy_number = max(min_copy, min(max_copy, copy_number))
    r["copy_number"] = copy_number

with out_tsv.open("w") as out:
    out.write("node\tlength\tdepth\tcopy_number\n")
    for r in rows:
        out.write(
            f'{r["node"]}\t{r["length"]}\t{r["depth"]:.10g}\t'
            f'{r["copy_number"]:.10g}\n'
        )

payload = {
    "mode": mode,
    "one_copy_depth": baseline,
    "fit_delta": fit_delta,
    "offset": offset,
    "min_copy": min_copy,
    "max_copy": max_copy,
    "cycle_nodes_excluded_from_initial_depth": sorted(cycle_nodes),
    "self_loop_nodes_excluded_from_initial_depth": sorted(self_loop_nodes),
}
with out_json.open("w") as out:
    json.dump(payload, out, indent=2)
