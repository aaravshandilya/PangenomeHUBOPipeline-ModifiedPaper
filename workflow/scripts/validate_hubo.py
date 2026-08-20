from pathlib import Path
import json
import itertools
import random

hubo_path = Path(snakemake.input.hubo)
out_path = Path(snakemake.output[0])
out_path.parent.mkdir(parents=True, exist_ok=True)

random_tests = int(snakemake.params.random_tests)
seed = int(snakemake.params.seed)
exhaustive_limit = int(snakemake.params.exhaustive_limit)
tol = float(snakemake.params.tolerance)

with hubo_path.open() as handle:
    data = json.load(handle)

variables = data["variables"]
terms = data["terms"]
meta = data["metadata"]

q = len(variables)
T = int(meta["walk_length_T"])
n_bits = int(meta["bits_per_time"])
N = int(meta["N_biological_nodes"])
lambda_edge = float(meta["lambda_edge"])
lambda_invalid = float(meta.get("lambda_invalid", 0.0))
weights = [float(x) for x in meta["node_weights"]]
edge_set = {tuple(e) for e in meta["oriented_edges"]}
state_labels = meta["state_labels"]

def polynomial_energy(bits):
    e = 0.0
    for term in terms:
        prod = 1
        for vid in term["variables"]:
            prod *= bits[vid]
            if prod == 0:
                break
        e += float(term["coefficient"]) * prod
    return e

def decode_states(bits):
    states = []
    for t in range(T):
        state = 0
        for k in range(n_bits):
            state |= (bits[t*n_bits+k] & 1) << k
        states.append(state)
    return states

def direct_energy(bits):
    states = decode_states(bits)
    h1 = 0.0
    for t in range(T - 1):
        if (states[t], states[t+1]) not in edge_set:
            h1 += lambda_edge

    hinvalid = sum(lambda_invalid for s in states if s >= 2*N)

    counts = [0] * N
    for state in states:
        if 0 <= state < 2*N:
            counts[state // 2] += 1
    h2 = sum((counts[i] - weights[i])**2 for i in range(N))
    return h1 + h2 + hinvalid

rng = random.Random(seed)
errors = []
examples = []
for _ in range(random_tests):
    bits = [rng.randint(0, 1) for _ in range(q)]
    ep = polynomial_energy(bits)
    ed = direct_energy(bits)
    err = abs(ep - ed)
    errors.append(err)
    if len(examples) < 5:
        examples.append({
            "polynomial_energy": ep,
            "direct_energy": ed,
            "absolute_error": err,
            "states": decode_states(bits),
        })

max_error = max(errors, default=0.0)
if max_error > tol:
    raise RuntimeError(
        f"HUBO validation failed: max direct-vs-polynomial error "
        f"{max_error} > tolerance {tol}"
    )

if q <= exhaustive_limit:
    best_energy = float("inf")
    best_bits = None
    for bits_tuple in itertools.product((0, 1), repeat=q):
        bits = list(bits_tuple)
        e = polynomial_energy(bits)
        if e < best_energy:
            best_energy = e
            best_bits = bits
    states = decode_states(best_bits)
    labels = [
        state_labels[s] if 0 <= s < len(state_labels) else f"INVALID_{s}"
        for s in states
    ]
    exhaustive = {
        "performed": True,
        "assignments_checked": 2 ** q,
        "best_energy": best_energy,
        "best_bits": best_bits,
        "best_states": states,
        "best_state_labels": labels,
    }
else:
    exhaustive = {
        "performed": False,
        "reason": f"q={q} exceeds exhaustive_limit={exhaustive_limit}",
    }

report = {
    "valid": True,
    "random_tests": random_tests,
    "maximum_absolute_error": max_error,
    "tolerance": tol,
    "binary_variable_count": q,
    "examples": examples,
    "exhaustive": exhaustive,
}
with out_path.open("w") as out:
    json.dump(report, out, indent=2)
