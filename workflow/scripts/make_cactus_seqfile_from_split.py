from pathlib import Path
import csv

split = Path(snakemake.input.split)
out = Path(snakemake.output[0])
out.parent.mkdir(parents=True, exist_ok=True)

rows=[]
with split.open(newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['role']=='train':
            rows.append(r)
if not rows:
    raise RuntimeError('No training genomes in split.tsv')

with out.open('w') as f:
    for r in rows:
        f.write(f"{r['sample']}\t{r['fasta']}\n")
