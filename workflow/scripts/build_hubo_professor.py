from pathlib import Path
import importlib.util
import csv, json, math

# Reuse the already validated Eq. (8)-(10) polynomial implementation.
# Snakemake executes script: files from a temporary location, so resolve the
# shared core from the workflow working directory instead of __file__.
base_path = Path.cwd() / "workflow" / "scripts" / "build_hubo.py"
if not base_path.exists():
    raise FileNotFoundError(f"Could not find HUBO core at {base_path}")
spec = importlib.util.spec_from_file_location('base_hubo', base_path)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

gfa = str(snakemake.input.gfa)
copies = str(snakemake.input.copies)

# Read weights to determine T using the policy requested in the study design.
weights=[]
with open(copies,newline='') as fh:
    for r in csv.DictReader(fh,delimiter='\t'):
        weights.append(float(r['copy_number']))
S=sum(weights)
policy=str(snakemake.params.fractional_T_policy)
if policy=='nearest': T=max(1,int(round(S)))
elif policy=='floor': T=max(1,int(math.floor(S)))
elif policy=='ceil': T=max(1,int(math.ceil(S)))
elif policy=='error':
    if abs(S-round(S))>1e-9: raise ValueError(f'sum(weights)={S} is not integral')
    T=max(1,int(round(S)))
else: raise ValueError(f'Unknown fractional_T_policy={policy}')

H, variables, meta = base.build_hubo(
    gfa,copies,
    lambda_edge=float(snakemake.params.lambda_edge),
    walk_length=T,
    max_variables=int(snakemake.params.max_variables),
    max_terms=int(snakemake.params.max_terms),
    tol=float(snakemake.params.coefficient_tolerance),
)

# Explicitly penalize unused binary state IDs X_t >= 2N.
# This is an extension requested for consistent treatment across graph builders.
N=int(meta['N_biological_nodes']); n_bits=int(meta['bits_per_time'])
lam_invalid=float(snakemake.params.lambda_invalid)
tol=float(snakemake.params.coefficient_tolerance)
max_terms=int(snakemake.params.max_terms)

def indicator(t,state):
    p=base.const_poly(1.0)
    for k in range(n_bits):
        x=base.var_poly(t*n_bits+k)
        b=(state>>k)&1
        factor=x if b else base.add_poly(base.const_poly(1.0),x,scale_b=-1.0,tol=tol)
        p=base.mul_poly(p,factor,tol=tol,max_terms=max_terms)
    return p

invalid_states=list(range(2*N,2**n_bits))
for t in range(T):
    for state in invalid_states:
        H=base.add_poly(H,indicator(t,state),scale_b=lam_invalid,tol=tol)
        if len(H)>max_terms: raise RuntimeError('HUBO exceeded max_terms while adding invalid-state penalty')

meta['lambda_invalid']=lam_invalid
meta['invalid_state_ids']=invalid_states
meta['walk_length_source']='sum_weights'
meta['sum_node_weights']=S
meta['fractional_T_policy']=policy
meta['term_count']=len(H)
meta['maximum_degree']=max((len(m) for m in H),default=0)

base.write_outputs(H,variables,meta,
    snakemake.output.hubo,snakemake.output.terms,
    snakemake.output.variable_map,snakemake.output.metadata)
