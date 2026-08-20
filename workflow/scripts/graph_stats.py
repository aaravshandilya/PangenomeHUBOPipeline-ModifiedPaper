from pathlib import Path
import json
from collections import Counter

gfa=Path(snakemake.input.gfa)
hubo=Path(snakemake.input.hubo)
out=Path(snakemake.output[0]); out.parent.mkdir(parents=True,exist_ok=True)

nodes=0; edges=0; lengths=[]; deg=Counter()
with gfa.open() as fh:
    for line in fh:
        f=line.rstrip('\n').split('\t')
        if not f: continue
        if f[0]=='S' and len(f)>=3:
            nodes+=1
            if f[2]!='*': lengths.append(len(f[2]))
            else:
                for t in f[3:]:
                    if t.startswith('LN:i:'): lengths.append(int(t.split(':')[-1])); break
        elif f[0]=='L' and len(f)>=5:
            edges+=1; deg[f[1]]+=1; deg[f[3]]+=1

with hubo.open() as fh: h=json.load(fh)
meta=h['metadata']; terms=h['terms']
degrees=[len(t['variables']) for t in terms]
coeffs=[float(t['coefficient']) for t in terms]
payload={
  'graph':{
    'nodes':nodes,'edges':edges,
    'mean_degree':(sum(deg.values())/nodes if nodes else 0),
    'max_degree':max(deg.values(),default=0),
    'mean_node_length':(sum(lengths)/len(lengths) if lengths else 0),
    'min_node_length':min(lengths,default=0),'max_node_length':max(lengths,default=0)
  },
  'hubo':{
    'N':meta['N_biological_nodes'],'T':meta['walk_length_T'],
    'bits_per_time':meta['bits_per_time'],'variables_q':meta['binary_variable_count'],
    'term_count':len(terms),'maximum_degree':max(degrees,default=0),
    'mean_degree':(sum(degrees)/len(degrees) if degrees else 0),
    'coefficient_min':min(coeffs,default=0),'coefficient_max':max(coeffs,default=0)
  }
}
with out.open('w') as fh: json.dump(payload,fh,indent=2)
