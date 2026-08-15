configfile: "config.yaml"

import csv
from pathlib import Path

SAMPLES_TSV = config["samples_tsv"]

with open(SAMPLES_TSV, newline="") as handle:
    SAMPLE_ROWS = list(csv.DictReader(handle, delimiter="\t"))

required_cols = {"sample", "fasta", "role"}
if not SAMPLE_ROWS or not required_cols.issubset(SAMPLE_ROWS[0].keys()):
    raise ValueError(
        f"{SAMPLES_TSV} must contain tab-separated columns: sample, fasta, role"
    )

HOLDOUT = config["holdout_sample"]
holdout_rows = [r for r in SAMPLE_ROWS if r["sample"] == HOLDOUT]
if len(holdout_rows) != 1:
    raise ValueError(f"Expected exactly one row for holdout_sample={HOLDOUT!r}")
HOLDOUT_FASTA = holdout_rows[0]["fasta"]

TRAINING_ROWS = [r for r in SAMPLE_ROWS if r["role"].lower() == "train"]
TRAINING_FASTAS = [r["fasta"] for r in TRAINING_ROWS]

REFERENCE = config["pangenome"]["reference"]
if REFERENCE not in {r["sample"] for r in TRAINING_ROWS}:
    raise ValueError(
        "pangenome.reference must name a training sample, not the held-out sample."
    )

PG_DIR = config["pangenome"]["out_dir"]
PG_NAME = config["pangenome"]["out_name"]
PG_GFA_GZ = f"{PG_DIR}/{PG_NAME}.full.gfa.gz"
PG_GBZ = f"{PG_DIR}/{PG_NAME}.full.gbz"

rule all:
    input:
        "results/hubo/hubo.json",
        "results/hubo/terms.tsv",
        "results/hubo/variable_map.tsv",
        "results/hubo/metadata.json",
        "results/hubo/validation.json",
        "results/subgraph/selected.gfa",
        "results/subgraph/selected_copy_numbers.tsv",
        "results/copy_number/copy_numbers.tsv",
        "results/annotation/annotated.gfa"

rule make_cactus_seqfile:
    input:
        manifest=SAMPLES_TSV,
        fastas=TRAINING_FASTAS
    output:
        "results/input/cactus.seqfile"
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/make_seqfile.py"

rule build_pangenome:
    input:
        seqfile="results/input/cactus.seqfile"
    output:
        gfa=PG_GFA_GZ,
        gbz=PG_GBZ
    params:
        outdir=PG_DIR,
        outname=PG_NAME,
        reference=REFERENCE,
        jobstore=config["pangenome"]["jobstore"],
        extra=config["pangenome"].get(
            "extra_args",
            "--noSplit --permissiveContigFilter"
        )
    threads:
        config["pangenome"].get("threads", 16)
    conda:
        "workflow/envs/cactus.yaml"
    log:
        "logs/minigraph_cactus.log"
    shell:
        r"""
        mkdir -p {params.outdir} "$(dirname {params.jobstore})"
        cactus-pangenome \
            {params.jobstore} \
            {input.seqfile} \
            --outDir {params.outdir} \
            --outName {params.outname} \
            --reference {params.reference} \
            --gfa full \
            --gbz full \
            --giraffe full \
            --mgCores {threads} \
            {params.extra} \
            > {log} 2>&1
        """

rule decompress_gfa:
    input:
        PG_GFA_GZ
    output:
        "results/pangenome/pangenome.full.gfa"
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        gzip -dc {input} > {output}
        """

rule simulate_reads:
    input:
        HOLDOUT_FASTA
    output:
        fastq=f"results/reads/{HOLDOUT}.fastq",
        truth=f"results/reads/{HOLDOUT}.truth.tsv"
    params:
        coverage=config["reads"]["coverage"],
        read_length=config["reads"]["read_length"],
        error_rate=config["reads"]["substitution_error_rate"],
        seed=config["reads"]["seed"]
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/simulate_reads.py"

rule map_reads:
    input:
        gbz=PG_GBZ,
        reads=f"results/reads/{HOLDOUT}.fastq"
    output:
        f"results/mapping/{HOLDOUT}.gaf"
    params:
        extra=config["mapping"].get("giraffe_args", "")
    threads:
        config["mapping"].get("threads", 8)
    conda:
        "workflow/envs/vg.yaml"
    log:
        "logs/giraffe.log"
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        vg giraffe \
            -t {threads} \
            -Z {input.gbz} \
            -f {input.reads} \
            -o gaf \
            {params.extra} \
            > {output} 2> {log}
        """

rule annotate_graph:
    input:
        gfa="results/pangenome/pangenome.full.gfa",
        gaf=f"results/mapping/{HOLDOUT}.gaf"
    output:
        gfa="results/annotation/annotated.gfa",
        depths="results/annotation/node_depths.tsv"
    params:
        min_mapq=config["annotation"]["min_mapq"]
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/annotate_gfa.py"

rule estimate_copy_numbers:
    input:
        gfa="results/annotation/annotated.gfa",
        depths="results/annotation/node_depths.tsv"
    output:
        copies="results/copy_number/copy_numbers.tsv",
        baseline="results/copy_number/baseline.json"
    params:
        mode=config["copy_number"]["mode"],
        min_depth=config["copy_number"]["min_depth"],
        depth_div=config["copy_number"]["depth_div"],
        offset=config["copy_number"]["offset"],
        min_copy=config["copy_number"]["min_copy"],
        max_copy=config["copy_number"]["max_copy"],
        fit_step=config["copy_number"]["fit_step"]
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/estimate_copy_numbers.py"

rule select_tangle:
    input:
        gfa="results/annotation/annotated.gfa",
        copies="results/copy_number/copy_numbers.tsv"
    output:
        gfa="results/subgraph/selected.gfa",
        copies="results/subgraph/selected_copy_numbers.tsv",
        summary="results/subgraph/summary.json"
    params:
        mode=config["subgraph"]["mode"],
        radius=config["subgraph"]["radius"],
        min_nodes=config["subgraph"]["min_nodes"],
        max_nodes=config["subgraph"]["max_nodes"],
        repeat_copy_threshold=config["subgraph"]["repeat_copy_threshold"]
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/select_tangle.py"

rule build_hubo:
    input:
        gfa="results/subgraph/selected.gfa",
        copies="results/subgraph/selected_copy_numbers.tsv"
    output:
        hubo="results/hubo/hubo.json",
        terms="results/hubo/terms.tsv",
        variable_map="results/hubo/variable_map.tsv",
        metadata="results/hubo/metadata.json"
    params:
        lambda_edge=config["hubo"]["lambda_edge"],
        walk_length=config["hubo"]["walk_length"],
        max_variables=config["hubo"]["max_variables"],
        max_terms=config["hubo"]["max_terms"],
        coefficient_tolerance=config["hubo"]["coefficient_tolerance"]
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/build_hubo.py"

rule validate_hubo:
    input:
        hubo="results/hubo/hubo.json"
    output:
        "results/hubo/validation.json"
    params:
        random_tests=config["validation"]["random_tests"],
        seed=config["validation"]["seed"],
        exhaustive_limit=config["validation"]["exhaustive_limit"],
        tolerance=config["validation"]["tolerance"]
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/validate_hubo.py"
