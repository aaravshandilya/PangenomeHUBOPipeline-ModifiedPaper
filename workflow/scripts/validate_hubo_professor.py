from pathlib import Path
import json,itertools,random

hubo_path=Path(snakemake.input.hubo); out_path=Path(snakemake.output[0]); out_path.parent.mkdir(parents=True,exist_ok=True)
with hubo_path.open() as fh: data=json.load(fh)
terms=data['terms']; meta=data['metadata']; variables=data['variables']
q=len(variables); T=int(meta['walk_length_T']); n_bits=int(meta['bits_per_time']); N=int(meta['N_biological_nodes'])
lambda_edge=float(meta['lambda_edge']); lambda_invalid=float(meta.get('lambda_invalid',0.0))
weights=[float(x) for x in meta['node_weights']]; edge_set={tuple(e) for e in meta['oriented_edges']}
labels=meta['state_labels']; tol=float(snakemake.params.tolerance)

def poly(bits):
    e=0.0
    for term in terms:
        p=1
        for v in term['variables']:
            p*=bits[v]
            if not p: break
        e+=float(term['coefficient'])*p
    return e

def decode(bits):
    out=[]
    for t in range(T):
        s=0
        for k in range(n_bits): s|=(bits[t*n_bits+k]&1)<<k
        out.append(s)
    return out

def direct(bits):
    st=decode(bits); h=0.0
    for s in st:
        if s>=2*N: h+=lambda_invalid
    for t in range(T-1):
        if (st[t],st[t+1]) not in edge_set: h+=lambda_edge
    counts=[0]*N
    for s in st:
        if 0<=s<2*N: counts[s//2]+=1
    h+=sum((counts[i]-weights[i])**2 for i in range(N))
    return h

rng=random.Random(int(snakemake.params.seed)); maxerr=0.0
for _ in range(int(snakemake.params.random_tests)):
    b=[rng.randint(0,1) for _ in range(q)]; maxerr=max(maxerr,abs(poly(b)-direct(b)))
if maxerr>tol: raise RuntimeError(f'HUBO validation failed: error {maxerr} > {tol}')

exhaustive={'performed':False,'reason':f'q={q} exceeds limit'}
limit=int(snakemake.params.exhaustive_limit)
if q<=limit:
    best=float('inf'); best_bits=None
    for bt in itertools.product((0,1),repeat=q):
        b=list(bt); e=poly(b)
        if e<best: best=e; best_bits=b
    st=decode(best_bits)
    exhaustive={'performed':True,'assignments_checked':2**q,'best_energy':best,'best_bits':best_bits,
                'best_states':st,'best_state_labels':[labels[s] if s<len(labels) else f'INVALID_{s}' for s in st]}

with out_path.open('w') as fh:
    json.dump({'valid':True,'maximum_absolute_error':maxerr,'binary_variable_count':q,'exhaustive':exhaustive},fh,indent=2)
