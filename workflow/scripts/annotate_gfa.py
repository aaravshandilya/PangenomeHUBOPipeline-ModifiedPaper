from collections import defaultdict
from pathlib import Path
import re

PATH_RE = re.compile(r"([<>])([^<>]+)")

def parse_tags(fields):
    tags = {}
    for x in fields:
        p = x.split(":", 2)
        if len(p) == 3:
            tags[p[0]] = (p[1], p[2])
    return tags

def invert_orient(o):
    return "-" if o == "+" else "+"

gfa_path = Path(snakemake.input.gfa)
gaf_path = Path(snakemake.input.gaf)
out_gfa = Path(snakemake.output.gfa)
out_depths = Path(snakemake.output.depths)
out_gfa.parent.mkdir(parents=True, exist_ok=True)

min_mapq = int(snakemake.params.min_mapq)

gfa_lines = []
lengths = {}
with gfa_path.open() as handle:
    for raw in handle:
        line = raw.rstrip("\n")
        gfa_lines.append(line)
        f = line.split("\t")
        if not f:
            continue
        if f[0] == "S":
            tags = parse_tags(f[3:])
            if len(f) >= 3 and f[2] != "*":
                lengths[f[1]] = len(f[2])
            elif "LN" in tags:
                lengths[f[1]] = int(float(tags["LN"][1]))
            else:
                raise ValueError(f"Segment {f[1]} has no sequence and no LN tag.")

# Keep one best primary-like alignment per read, reproducing the intent of
# qpg/tag_gfa_ga.pl while being robust to non-adjacent GAF records.
alignments = defaultdict(list)
with gaf_path.open() as handle:
    for raw in handle:
        if not raw.strip() or raw.startswith("#"):
            continue
        f = raw.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        try:
            mapq = int(f[11])
            block_len = int(f[10])
        except ValueError:
            continue
        if mapq < min_mapq:
            continue
        tags = parse_tags(f[12:])
        tp = tags.get("tp", (None, None))[1]
        is_secondary = (tp == "S")
        alignments[f[0]].append((is_secondary, block_len, mapq, f))

coverage = defaultdict(float)
edge_counts = defaultdict(int)

for qname, candidates in alignments.items():
    primaries = [x for x in candidates if not x[0]]
    pool = primaries if primaries else candidates
    _, _, _, f = max(pool, key=lambda x: (x[1], x[2]))

    path_field = f[5]
    path = PATH_RE.findall(path_field)
    if not path:
        continue

    # GAF query coordinates are 0-based half-open.
    try:
        qstart, qend = int(f[2]), int(f[3])
        aligned_query_len = max(0, qend - qstart)
    except ValueError:
        aligned_query_len = int(f[10])

    nodes = [node for _, node in path]
    orients = ["+" if marker == ">" else "-" for marker, _ in path]

    # Edge traversal counts.
    for a in range(len(nodes) - 1):
        edge_counts[(nodes[a], orients[a], nodes[a+1], orients[a+1])] += 1

    if len(nodes) == 1:
        n = nodes[0]
        coverage[n] += min(lengths.get(n, aligned_query_len), aligned_query_len)
        continue

    # qpg/tag_gfa_ga.pl gives internal graph nodes their full segment length
    # and distributes the remainder between the first and last nodes.
    internal_used = 0.0
    for n in nodes[1:-1]:
        ln = lengths[n]
        coverage[n] += ln
        internal_used += ln

    remaining = max(0.0, aligned_query_len - internal_used)
    first = nodes[0]
    last = nodes[-1]

    first_used = max(1.0, remaining / 2.0) if remaining > 0 else 0.0
    first_used = min(first_used, lengths[first])
    coverage[first] += first_used

    last_used = max(0.0, remaining - first_used)
    last_used = min(last_used, lengths[last])
    coverage[last] += last_used

def replace_tags(fields, replacements):
    remove = set(replacements)
    kept = []
    for x in fields:
        key = x.split(":", 1)[0]
        if key not in remove:
            kept.append(x)
    return kept + list(replacements.values())

with out_gfa.open("w") as out:
    for line in gfa_lines:
        f = line.split("\t")
        if not f:
            continue
        if f[0] == "S":
            node = f[1]
            ln = lengths[node]
            cov_bases = coverage.get(node, 0.0)
            depth = cov_bases / ln if ln else 0.0
            prefix = f[:3]
            tags = replace_tags(
                f[3:],
                {
                    "LN": f"LN:i:{ln}",
                    "KC": f"KC:i:{int(round(cov_bases))}",
                    "SC": f"SC:f:{depth:.10g}",
                },
            )
            out.write("\t".join(prefix + tags) + "\n")
        elif f[0] == "L" and len(f) >= 6:
            u, ou, v, ov = f[1], f[2], f[3], f[4]
            direct = edge_counts.get((u, ou, v, ov), 0)
            reverse = edge_counts.get(
                (v, invert_orient(ov), u, invert_orient(ou)), 0
            )
            ec = direct if (u, ou, v, ov) == (
                v, invert_orient(ov), u, invert_orient(ou)
            ) else direct + reverse
            tags = replace_tags(f[6:], {"EC": f"EC:i:{ec}"})
            out.write("\t".join(f[:6] + tags) + "\n")
        else:
            out.write(line + "\n")

with out_depths.open("w") as out:
    out.write("node\tlength\tcovered_bases\tdepth\n")
    for node, ln in lengths.items():
        cov_bases = coverage.get(node, 0.0)
        out.write(
            f"{node}\t{ln}\t{cov_bases:.8f}\t"
            f"{(cov_bases / ln if ln else 0.0):.10g}\n"
        )
