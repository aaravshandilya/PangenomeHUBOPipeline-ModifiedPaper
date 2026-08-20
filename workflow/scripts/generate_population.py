from pathlib import Path
import shlex, subprocess

exe = str(snakemake.input.exe)
out_fa = Path(snakemake.output.population_fasta)
out_manifest = Path(snakemake.output.manifest)
out_fa.parent.mkdir(parents=True, exist_ok=True)

cmd = [exe] + shlex.split(str(snakemake.params.genome_opts)) + ["-s", str(snakemake.params.seed)]
with out_fa.open("w") as out:
    subprocess.run(cmd, stdout=out, check=True)

records = []
name = None
seq = []
with out_fa.open() as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
if name is not None:
    records.append((name, "".join(seq)))

if len(records) < 50:
    raise RuntimeError(f"Expected at least 50 genomes, found {len(records)}")

individual_dir = out_fa.parent / "genomes"
individual_dir.mkdir(exist_ok=True)
with out_manifest.open("w") as mf:
    mf.write("index\tsample\tfasta\tlength\n")
    for i, (name, sequence) in enumerate(records):
        path = individual_dir / f"{name}.fa"
        with path.open("w") as f:
            f.write(f">{name}\n")
            for j in range(0, len(sequence), 80):
                f.write(sequence[j:j+80] + "\n")
        mf.write(f"{i}\t{name}\t{path.resolve()}\t{len(sequence)}\n")

print(f"Generated {len(records)} synthetic genomes with seed {snakemake.params.seed}")
