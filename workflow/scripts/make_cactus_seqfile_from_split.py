from pathlib import Path
import csv
import re
import shutil

split = Path(snakemake.input.split)
out = Path(snakemake.output[0])
out.parent.mkdir(parents=True, exist_ok=True)

rows = []
with split.open(newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["role"] == "train":
            rows.append(r)

if not rows:
    raise RuntimeError("No training genomes in split.tsv")

# Cactus seqfiles treat '#' as the start of a comment. The qpg-generated
# sample/file names contain '#', so stage byte-identical FASTAs under safe
# names before writing the seqfile. This changes only labels/paths, not the
# biological sequences supplied to the graph builder.
staging = out.parent / "input_fastas"
if staging.exists():
    shutil.rmtree(staging)
staging.mkdir(parents=True, exist_ok=True)

used = set()
with out.open("w") as f:
    for i, r in enumerate(rows):
        safe_sample = re.sub(r"[^A-Za-z0-9_.-]+", "_", r["sample"]).strip("_")
        if not safe_sample:
            safe_sample = f"sample_{i:03d}"
        base = safe_sample
        j = 1
        while safe_sample in used:
            safe_sample = f"{base}_{j}"
            j += 1
        used.add(safe_sample)

        src = Path(r["fasta"]).resolve()
        if not src.exists():
            raise FileNotFoundError(f"Training FASTA not found: {src}")

        dst = staging / f"{safe_sample}.fa"
        shutil.copyfile(src, dst)
        f.write(f"{safe_sample}\t{dst.resolve()}\n")
