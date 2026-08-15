import csv
from pathlib import Path

manifest = Path(snakemake.input.manifest)
out = Path(snakemake.output[0])
out.parent.mkdir(parents=True, exist_ok=True)

rows = []
with manifest.open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["role"].strip().lower() == "train":
            fasta = Path(row["fasta"]).resolve()
            if not fasta.exists():
                raise FileNotFoundError(f"Training FASTA not found: {fasta}")
            rows.append((row["sample"].strip(), fasta))

if len(rows) < 2:
    raise ValueError("Use at least two training assemblies for a pangenome.")

with out.open("w") as handle:
    for sample, fasta in rows:
        handle.write(f"{sample}\t{fasta}\n")
