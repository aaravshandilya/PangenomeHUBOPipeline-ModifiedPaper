"""
Build the oriented-tangle HUBO from Cudby & Strelchuk (2026), Eqs. (8)-(10).

The binary polynomial is represented as:
    monomial frozenset({variable_ids}) -> coefficient

Because all variables are binary, x_i^2 = x_i, so multiplying monomials
corresponds to taking the union of their variable sets.
"""

from collections import defaultdict
from pathlib import Path
import argparse
import csv
import json
import math

def clean(poly, tol=1e-12):
    return {m: c for m, c in poly.items() if abs(c) > tol}

def add_poly(a, b, scale_b=1.0, tol=1e-12):
    out = defaultdict(float)
    out.update(a)
    for m, c in b.items():
        out[m] += scale_b * c
    return clean(out, tol)

def scale_poly(a, scale, tol=1e-12):
    return clean({m: scale * c for m, c in a.items()}, tol)

def mul_poly(a, b, tol=1e-12, max_terms=None):
    out = defaultdict(float)
    for ma, ca in a.items():
        for mb, cb in b.items():
            out[ma | mb] += ca * cb
    out = clean(out, tol)
    if max_terms is not None and len(out) > max_terms:
        raise RuntimeError(
            f"Polynomial expanded to {len(out):,} terms, exceeding "
            f"max_terms={max_terms:,}. Select a smaller tangle or lower T."
        )
    return out

def const_poly(c):
    return {frozenset(): float(c)}

def var_poly(var_id):
    return {frozenset({int(var_id)}): 1.0}

def parse_gfa(path):
    node_order = []
    node_lengths = {}
    raw_links = []
    with open(path) as handle:
        for raw in handle:
            f = raw.rstrip("\n").split("\t")
            if not f:
                continue
            if f[0] == "S":
                node_order.append(f[1])
                if len(f) >= 3 and f[2] != "*":
                    node_lengths[f[1]] = len(f[2])
                else:
                    ln = None
                    for tag in f[3:]:
                        if tag.startswith("LN:i:"):
                            ln = int(tag.split(":")[-1])
                    if ln is None:
                        raise ValueError(f"No sequence/LN for node {f[1]}")
                    node_lengths[f[1]] = ln
            elif f[0] == "L" and len(f) >= 5:
                raw_links.append((f[1], f[2], f[3], f[4]))
    return node_order, node_lengths, raw_links

def load_weights(path):
    weights = {}
    with open(path, newline="") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            weights[r["node"]] = float(r["copy_number"])
    return weights

def state_id(node_index, orient):
    return 2 * node_index + (0 if orient == "+" else 1)

def inv_orient(o):
    return "-" if o == "+" else "+"

def build_hubo(
    gfa_path,
    copy_path,
    lambda_edge=10.0,
    walk_length="auto",
    max_variables=40,
    max_terms=2_000_000,
    tol=1e-12,
):
    node_order, node_lengths, links = parse_gfa(gfa_path)
    weights_by_name = load_weights(copy_path)

    if not node_order:
        raise ValueError("Selected subgraph contains no nodes.")
    missing = [n for n in node_order if n not in weights_by_name]
    if missing:
        raise ValueError(f"Missing copy numbers for nodes: {missing}")

    N = len(node_order)
    n_bits = math.ceil(math.log2(2 * N))
    n_bits = max(1, n_bits)

    biological_weights = [weights_by_name[n] for n in node_order]
    if str(walk_length).lower() == "auto":
        T = max(1, int(round(sum(biological_weights))))
    else:
        T = max(1, int(walk_length))

    q = T * n_bits
    if q > int(max_variables):
        raise RuntimeError(
            f"HUBO would require q=T*ceil(log2(2N))={T}*{n_bits}={q} "
            f"binary variables, exceeding max_variables={max_variables}. "
            "Reduce the selected subgraph/walk length or raise the guard deliberately."
        )

    node_index = {name: i for i, name in enumerate(node_order)}

    # Directed oriented edge set. GFA links imply a reverse-complement edge too.
    oriented_edges = set()
    for u, ou, v, ov in links:
        if u not in node_index or v not in node_index:
            continue
        iu = state_id(node_index[u], ou)
        iv = state_id(node_index[v], ov)
        oriented_edges.add((iu, iv))
        oriented_edges.add((
            state_id(node_index[v], inv_orient(ov)),
            state_id(node_index[u], inv_orient(ou)),
        ))

    def var_id(t, k):
        return t * n_bits + k

    indicator_cache = {}

    def indicator(t, state):
        """
        Eq. (8):
          1(X_t=i) = prod_k (1-b_k-x_tk+2 b_k x_tk)

        For b_k=0 the factor is (1-x_tk); for b_k=1 it is x_tk.
        """
        key = (t, state)
        if key in indicator_cache:
            return indicator_cache[key]
        p = const_poly(1.0)
        for k in range(n_bits):
            b = (state >> k) & 1
            x = var_poly(var_id(t, k))
            factor = x if b else add_poly(const_poly(1.0), x, scale_b=-1.0, tol=tol)
            p = mul_poly(p, factor, tol=tol, max_terms=max_terms)
        indicator_cache[key] = p
        return p

    H = {}

    # Eq. (9): edge-following penalty.
    # H1 = Lambda1 * sum_t [1 - sum_(i,j in E) I(X_t=i) I(X_t+1=j)]
    for t in range(T - 1):
        H = add_poly(H, const_poly(lambda_edge), tol=tol)
        for i, j in oriented_edges:
            product_ij = mul_poly(
                indicator(t, i),
                indicator(t + 1, j),
                tol=tol,
                max_terms=max_terms,
            )
            H = add_poly(H, product_ij, scale_b=-lambda_edge, tol=tol)
            if len(H) > max_terms:
                raise RuntimeError(
                    f"HUBO exceeded max_terms={max_terms:,} while building H1."
                )

    # Eq. (10): copy-number / visitation objective.
    # For biological node l, positive and negative states are 2l and 2l+1.
    for l, weight in enumerate(biological_weights):
        s = const_poly(-weight)
        for t in range(T):
            s = add_poly(s, indicator(t, 2*l), tol=tol)
            s = add_poly(s, indicator(t, 2*l + 1), tol=tol)
        squared = mul_poly(s, s, tol=tol, max_terms=max_terms)
        H = add_poly(H, squared, tol=tol)
        if len(H) > max_terms:
            raise RuntimeError(
                f"HUBO exceeded max_terms={max_terms:,} while building H2."
            )

    H = clean(H, tol)

    variables = []
    for t in range(T):
        for k in range(n_bits):
            vid = var_id(t, k)
            variables.append({
                "id": vid,
                "name": f"x_{t}_{k}",
                "time": t,
                "bit": k,
            })

    state_labels = []
    for l, node in enumerate(node_order):
        state_labels.extend([f"{node}+", f"{node}-"])

    max_degree = max((len(m) for m in H), default=0)
    degree_hist = defaultdict(int)
    for m in H:
        degree_hist[len(m)] += 1

    metadata = {
        "formulation": "oriented_tangle_hubo_cudby_strelchuk_2026",
        "N_biological_nodes": N,
        "oriented_state_count": 2 * N,
        "bits_per_time": n_bits,
        "walk_length_T": T,
        "binary_variable_count": q,
        "lambda_edge": float(lambda_edge),
        "node_order": node_order,
        "node_lengths": node_lengths,
        "node_weights": biological_weights,
        "state_labels": state_labels,
        "oriented_edges": [list(e) for e in sorted(oriented_edges)],
        "term_count": len(H),
        "maximum_degree": max_degree,
        "degree_histogram": {str(k): v for k, v in sorted(degree_hist.items())},
        "coefficient_tolerance": tol,
    }
    return H, variables, metadata

def write_outputs(H, variables, metadata, hubo_path, terms_path, varmap_path, meta_path):
    id_to_name = {v["id"]: v["name"] for v in variables}
    terms = [
        {
            "variables": sorted(m),
            "coefficient": c,
        }
        for m, c in sorted(
            H.items(),
            key=lambda item: (len(item[0]), tuple(sorted(item[0])))
        )
    ]

    payload = {
        "format": "binary_polynomial_v1",
        "variables": variables,
        "terms": terms,
        "metadata": metadata,
    }

    Path(hubo_path).parent.mkdir(parents=True, exist_ok=True)
    with open(hubo_path, "w") as out:
        json.dump(payload, out, indent=2)

    with open(meta_path, "w") as out:
        json.dump(metadata, out, indent=2)

    with open(varmap_path, "w") as out:
        out.write("id\tname\ttime\tbit\n")
        for v in variables:
            out.write(f'{v["id"]}\t{v["name"]}\t{v["time"]}\t{v["bit"]}\n')

    with open(terms_path, "w") as out:
        out.write("degree\tcoefficient\tvariables\n")
        for term in terms:
            names = ",".join(id_to_name[i] for i in term["variables"])
            out.write(
                f'{len(term["variables"])}\t{term["coefficient"]:.17g}\t{names}\n'
            )

def run_from_snakemake(sm):
    H, variables, metadata = build_hubo(
        sm.input.gfa,
        sm.input.copies,
        lambda_edge=float(sm.params.lambda_edge),
        walk_length=sm.params.walk_length,
        max_variables=int(sm.params.max_variables),
        max_terms=int(sm.params.max_terms),
        tol=float(sm.params.coefficient_tolerance),
    )
    write_outputs(
        H, variables, metadata,
        sm.output.hubo, sm.output.terms, sm.output.variable_map, sm.output.metadata
    )

def cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--hubo", required=True)
    ap.add_argument("--terms", required=True)
    ap.add_argument("--variable-map", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--lambda-edge", type=float, default=10.0)
    ap.add_argument("--walk-length", default="auto")
    ap.add_argument("--max-variables", type=int, default=40)
    ap.add_argument("--max-terms", type=int, default=2_000_000)
    ap.add_argument("--tol", type=float, default=1e-12)
    args = ap.parse_args()
    H, variables, metadata = build_hubo(
        args.gfa, args.copies, args.lambda_edge, args.walk_length,
        args.max_variables, args.max_terms, args.tol
    )
    write_outputs(
        H, variables, metadata,
        args.hubo, args.terms, args.variable_map, args.metadata
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
elif __name__ == "__main__":
    cli()
