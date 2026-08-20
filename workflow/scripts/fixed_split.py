from pathlib import Path
import random

pop = Path(snakemake.input.pop)
out_split = Path(snakemake.output.split)
out_train = Path(snakemake.output.train_fofn)
out_test = Path(snakemake.output.test_fofn)
out_heldout = Path(snakemake.output.heldout)
out_concat = Path(snakemake.output.train_concat)
for p in [out_split, out_train, out_test, out_heldout, out_concat]:
    p.parent.mkdir(parents=True, exist_ok=True)

records=[]; name=None; seq=[]
with pop.open() as fh:
    for line in fh:
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if name is not None: records.append((name,''.join(seq)))
            name=line[1:].split()[0]; seq=[]
        else:
            seq.append(line)
if name is not None: records.append((name,''.join(seq)))

# Match the original qpg study design closely: choose from the last 50 genomes,
# but make the shuffle reproducible with a documented seed.
candidates = records[-50:]
rng = random.Random(int(snakemake.params.split_seed))
rng.shuffle(candidates)
n_train = int(snakemake.params.n_training)
n_test = int(snakemake.params.n_test)
if n_train + n_test > len(candidates):
    raise ValueError('Requested split exceeds 50-genome candidate pool')
train = candidates[:n_train]
test = candidates[n_train:n_train+n_test]
hidx = int(snakemake.params.heldout_index)
heldout = test[hidx]

gdir = pop.parent/'genomes'
def path_for(name): return (gdir/f'{name}.fa').resolve()

with out_train.open('w') as f:
    for n,_ in train: f.write(str(path_for(n))+'\n')
with out_test.open('w') as f:
    for n,_ in test: f.write(str(path_for(n))+'\n')
with out_split.open('w') as f:
    f.write('sample\trole\tfasta\n')
    for n,_ in train: f.write(f'{n}\ttrain\t{path_for(n)}\n')
    for n,_ in test: f.write(f'{n}\ttest\t{path_for(n)}\n')

with out_heldout.open('w') as f:
    f.write(f'>{heldout[0]}\n')
    s=heldout[1]
    for i in range(0,len(s),80): f.write(s[i:i+80]+'\n')

# PGGB requires globally unique path names. Prefix each sequence header with sample.
with out_concat.open('w') as out:
    for n,s in train:
        out.write(f'>{n}#1#contig\n')
        for i in range(0,len(s),80): out.write(s[i:i+80]+'\n')

print(f'Training genomes: {len(train)}; test genomes: {len(test)}; held out: {heldout[0]}')
