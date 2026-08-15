from pathlib import Path
import math
import random

DNA = "ACGT"

def read_fasta(path):
    name = None
    chunks = []
    records = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return records

def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]

def mutate_substitutions(seq, error_rate, rng):
    out = []
    for base in seq:
        b = base if base in DNA else "N"
        if b in DNA and rng.random() < error_rate:
            choices = [x for x in DNA if x != b]
            b = rng.choice(choices)
        out.append(b)
    return "".join(out)

fasta = Path(snakemake.input[0])
fastq = Path(snakemake.output.fastq)
truth = Path(snakemake.output.truth)
fastq.parent.mkdir(parents=True, exist_ok=True)

coverage = float(snakemake.params.coverage)
read_len = int(snakemake.params.read_length)
error_rate = float(snakemake.params.error_rate)
rng = random.Random(int(snakemake.params.seed))

records = [(n, s) for n, s in read_fasta(fasta) if len(s) >= read_len]
if not records:
    raise ValueError(
        f"No FASTA record in {fasta} is at least read_length={read_len} bp."
    )

total_bases = sum(len(s) for _, s in records)
n_reads = max(1, math.ceil(coverage * total_bases / read_len))

weights = [len(s) for _, s in records]
with fastq.open("w") as fq, truth.open("w") as tr:
    tr.write("read\tcontig\tstart\tend\tstrand\n")
    for idx in range(n_reads):
        name, seq = rng.choices(records, weights=weights, k=1)[0]
        start = rng.randrange(0, len(seq) - read_len + 1)
        raw = seq[start:start + read_len]
        strand = "+" if rng.random() < 0.5 else "-"
        observed = raw if strand == "+" else revcomp(raw)
        observed = mutate_substitutions(observed, error_rate, rng)
        qname = f"read_{idx:08d}"
        fq.write(f"@{qname}\n{observed}\n+\n{'I' * len(observed)}\n")
        tr.write(f"{qname}\t{name}\t{start}\t{start + read_len}\t{strand}\n")
