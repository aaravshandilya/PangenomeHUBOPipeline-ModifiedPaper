configfile: "config.yaml"

from pathlib import Path

BUILDERS = config["builders"]["names"]
SEED = int(config["qpg"]["population_seed"])

RAW_GFA = {
    "minigraph": "results/graphs/minigraph/raw.gfa",
    "minigraph_cactus": "results/graphs/minigraph_cactus/raw.gfa",
    "pggb": "results/graphs/pggb/raw.gfa",
}

rule all:
    input:
        expand("results/{builder}/hubo/hubo.json", builder=BUILDERS),
        expand("results/{builder}/hubo/hubo_terms.tsv", builder=BUILDERS),
        expand("results/{builder}/hubo/metadata.json", builder=BUILDERS),
        expand("results/{builder}/hubo/validation.json", builder=BUILDERS),
        expand("results/{builder}/graph_stats.json", builder=BUILDERS),
        "results/population/split.tsv",
        "results/reads/heldout.fastq"

rule prepare_qpg:
    output:
        exe="resources/qpg/genome_create",
        commit="resources/qpg/COMMIT"
    params:
        repo=config["qpg"]["repo_url"],
        commit=config["qpg"]["commit"]
    conda:
        "workflow/envs/build.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p resources
        if [ ! -d resources/qpg/.git ]; then
            git clone {params.repo} resources/qpg
        fi
        git -C resources/qpg fetch --all --tags
        git -C resources/qpg checkout {params.commit}
        gcc -O2 resources/qpg/genome_create.c -o {output.exe} -lm
        git -C resources/qpg rev-parse HEAD > {output.commit}
        """

rule generate_population:
    input:
        exe=rules.prepare_qpg.output.exe
    output:
        population_fasta="results/population/population.fa",
        manifest="results/population/population.tsv"
    params:
        genome_opts=config["qpg"]["genome_opts"],
        seed=SEED
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/generate_population.py"

rule fixed_split:
    input:
        population_fasta=rules.generate_population.output.population_fasta
    output:
        split="results/population/split.tsv",
        train_fofn="results/population/train.fofn",
        test_fofn="results/population/test.fofn",
        heldout="results/population/heldout.fa",
        train_concat="results/population/training.fa"
    params:
        n_training=config["population"]["n_training"],
        n_test=config["population"]["n_test"],
        split_seed=config["population"]["split_seed"],
        heldout_index=config["population"]["heldout_index"]
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/fixed_split.py"

rule cactus_seqfile:
    input:
        split=rules.fixed_split.output.split
    output:
        "results/graphs/minigraph_cactus/seqfile.txt"
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/make_cactus_seqfile_from_split.py"

rule build_minigraph:
    input:
        fofn=rules.fixed_split.output.train_fofn
    output:
        RAW_GFA["minigraph"]
    threads:
        config["builders"]["minigraph"]["threads"]
    params:
        args=config["builders"]["minigraph"]["args"]
    conda:
        "workflow/envs/minigraph.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output})
        minigraph {params.args} $(cat {input.fofn}) > {output}
        """

rule build_minigraph_cactus:
    input:
        seqfile=rules.cactus_seqfile.output
    output:
        RAW_GFA["minigraph_cactus"]
    threads:
        config["builders"]["minigraph_cactus"]["threads"]
    params:
        args=config["builders"]["minigraph_cactus"]["args"]
    conda:
        "workflow/envs/cactus.yaml"
    shell:
        r"""
        set -euo pipefail
        outdir=results/graphs/minigraph_cactus/cactus_out
        jobstore=results/graphs/minigraph_cactus/jobstore
        rm -rf "$jobstore" "$outdir"
        mkdir -p "$outdir"
        ref=$(head -1 {input.seqfile} | cut -f1)
        cactus-pangenome "$jobstore" {input.seqfile} \
            --outDir "$outdir" --outName graph --reference "$ref" \
            --gfa full --mgCores {threads} {params.args}
        gfa=$(find "$outdir" -type f \( -name '*.full.gfa' -o -name '*.full.gfa.gz' \) | head -1)
        test -n "$gfa"
        if [[ "$gfa" == *.gz ]]; then gzip -dc "$gfa" > {output}; else cp "$gfa" {output}; fi
        """

rule build_pggb:
    input:
        fasta=rules.fixed_split.output.train_concat,
        fofn=rules.fixed_split.output.train_fofn
    output:
        RAW_GFA["pggb"]
    threads:
        config["builders"]["pggb"]["threads"]
    params:
        p=config["builders"]["pggb"]["map_pct_id"],
        s=config["builders"]["pggb"]["segment_length"]
    conda:
        "workflow/envs/pggb.yaml"
    shell:
        r"""
        set -euo pipefail
        outdir=results/graphs/pggb/pggb_out
        rm -rf "$outdir" && mkdir -p "$outdir"
        n=$(wc -l < {input.fofn})
        pggb -i {input.fasta} -o "$outdir" -n "$n" -t {threads} -p {params.p} -s {params.s}
        gfa=$(find "$outdir" -type f -name '*.gfa' | head -1)
        test -n "$gfa"
        cp "$gfa" {output}
        """

rule normalize_graph:
    input:
        lambda wc: RAW_GFA[wc.builder]
    output:
        gfa="results/{builder}/normalized/graph.gfa",
        nodes="results/{builder}/normalized/nodes.tsv",
        edges="results/{builder}/normalized/edges.tsv",
        seqs="results/{builder}/normalized/node_sequences.fa"
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/normalize_graph.py"

rule simulate_reads_once:
    input:
        rules.fixed_split.output.heldout
    output:
        fastq="results/reads/heldout.fastq",
        truth="results/reads/heldout.truth.tsv"
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
        gfa="results/{builder}/normalized/graph.gfa",
        reads=rules.simulate_reads_once.output.fastq
    output:
        "results/{builder}/mapping/reads.gaf"
    threads:
        config["mapping"]["threads"]
    params:
        args=config["mapping"]["graphaligner_args"]
    conda:
        "workflow/envs/graphaligner.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output})
        GraphAligner -g {input.gfa} -f {input.reads} -a {output} -t {threads} {params.args}
        """

rule annotate_graph:
    input:
        gfa="results/{builder}/normalized/graph.gfa",
        gaf="results/{builder}/mapping/reads.gaf"
    output:
        gfa="results/{builder}/weighted/annotated.gfa",
        depths="results/{builder}/weighted/node_depths.tsv"
    params:
        min_mapq=config["mapping"]["min_mapq"]
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/annotate_gfa.py"

rule estimate_copy_numbers:
    input:
        gfa="results/{builder}/weighted/annotated.gfa",
        depths="results/{builder}/weighted/node_depths.tsv"
    output:
        copies="results/{builder}/weighted/weights.tsv",
        baseline="results/{builder}/weighted/baseline.json"
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
        gfa="results/{builder}/weighted/annotated.gfa",
        copies="results/{builder}/weighted/weights.tsv"
    output:
        gfa="results/{builder}/tangle/selected.gfa",
        copies="results/{builder}/tangle/weights.tsv",
        summary="results/{builder}/tangle/summary.json"
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
        gfa="results/{builder}/tangle/selected.gfa",
        copies="results/{builder}/tangle/weights.tsv"
    output:
        hubo="results/{builder}/hubo/hubo.json",
        terms="results/{builder}/hubo/hubo_terms.tsv",
        variable_map="results/{builder}/hubo/variable_map.tsv",
        metadata="results/{builder}/hubo/metadata.json"
    params:
        lambda_edge=config["hubo"]["lambda_edge"],
        lambda_invalid=config["hubo"]["lambda_invalid"],
        walk_length=config["hubo"]["walk_length"],
        fractional_T_policy=config["hubo"]["fractional_T_policy"],
        max_variables=config["hubo"]["max_variables"],
        max_terms=config["hubo"]["max_terms"],
        coefficient_tolerance=config["hubo"]["coefficient_tolerance"]
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/build_hubo_professor.py"

rule validate_hubo:
    input:
        hubo="results/{builder}/hubo/hubo.json"
    output:
        "results/{builder}/hubo/validation.json"
    params:
        random_tests=config["validation"]["random_tests"],
        seed=config["validation"]["seed"],
        exhaustive_limit=config["validation"]["exhaustive_limit"],
        tolerance=config["validation"]["tolerance"]
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/validate_hubo.py"

rule graph_stats:
    input:
        gfa="results/{builder}/normalized/graph.gfa",
        hubo="results/{builder}/hubo/hubo.json"
    output:
        "results/{builder}/graph_stats.json"
    conda:
        "workflow/envs/python.yaml"
    script:
        "workflow/scripts/graph_stats.py"
