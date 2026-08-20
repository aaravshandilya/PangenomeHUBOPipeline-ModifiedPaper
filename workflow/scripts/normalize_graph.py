from pathlib import Path

src = Path(snakemake.input[0])
out_gfa = Path(snakemake.output.gfa)
out_nodes = Path(snakemake.output.nodes)
out_edges = Path(snakemake.output.edges)
out_seqs = Path(snakemake.output.seqs)
for p in [out_gfa,out_nodes,out_edges,out_seqs]: p.parent.mkdir(parents=True, exist_ok=True)

nodes=[]; links=[]; headers=[]
with src.open() as fh:
    for raw in fh:
        line=raw.rstrip('\n')
        if not line: continue
        f=line.split('\t')
        if f[0]=='H': headers.append(line)
        elif f[0]=='S' and len(f)>=3:
            nodes.append((f[1],f[2],f[3:]))
        elif f[0]=='L' and len(f)>=6:
            links.append((f[1],f[2],f[3],f[4],f[5],f[6:]))

if not nodes:
    raise RuntimeError(f'No GFA segment records found in {src}')

with out_gfa.open('w') as out:
    out.write((headers[0] if headers else 'H\tVN:Z:1.0')+'\n')
    for n,seq,tags in nodes:
        out.write('\t'.join(['S',n,seq]+tags)+'\n')
    for u,ou,v,ov,ovlp,tags in links:
        out.write('\t'.join(['L',u,ou,v,ov,ovlp]+tags)+'\n')

with out_nodes.open('w') as out:
    out.write('node\tlength\n')
    for n,seq,tags in nodes:
        ln = len(seq) if seq!='*' else next((int(t.split(':')[-1]) for t in tags if t.startswith('LN:i:')),0)
        out.write(f'{n}\t{ln}\n')

with out_edges.open('w') as out:
    out.write('from\tfrom_orient\tto\tto_orient\n')
    for u,ou,v,ov,_,_ in links:
        out.write(f'{u}\t{ou}\t{v}\t{ov}\n')

with out_seqs.open('w') as out:
    for n,seq,_ in nodes:
        if seq=='*': continue
        out.write(f'>{n}\n')
        for i in range(0,len(seq),80): out.write(seq[i:i+80]+'\n')

print(f'Normalized {len(nodes)} nodes and {len(links)} edges from {src}')
