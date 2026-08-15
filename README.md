# Snakemake Minigraph-Cactus → PGSA → HUBO pipeline

This workflow recreates the **structure** of the PGSA pipeline from Cudby et al.
(2026), but replaces its QUBO stage with the **oriented-tangle HUBO** introduced
by Cudby & Strelchuk (2026).

## What the workflow does

```text
training assemblies
       |
       v
Minigraph-Cactus
       |
       +---- full GFA
       |
       v
current vg autoindex (sr-giraffe)
       |
       +---- Giraffe GBZ + mapping indexes
       |
held-out assembly
       |
       v
simulate ~30x single-end short reads
       |
       v
vg Giraffe -> GAF
       |
       v
Python graph annotation
  LN / KC / SC / EC
       |
       v
Python copy-number estimation w(v)
       |
       v
select small branch/repeat-rich tangle
       |
       v
build oriented HUBO in Python
       |
       v
validate polynomial against direct objective
```

The workflow intentionally **stops after generating and validating the HUBO**.
A solver can be added later without changing the biological preprocessing.

## Relationship to the original PGSA code

The original `jkbonfield/qpg` workflow uses shell/Perl/Python scripts. This
Snakemake implementation preserves the major steps but makes the custom stages
Python rules:

| PGSA concept | This workflow |
|---|---|
| Create pangenome | Minigraph-Cactus |
| Hold out target genome | `samples.tsv` role=`test` |
| Simulate shotgun reads | `simulate_reads.py` |
| Map reads to graph | `vg giraffe` |
| Add node/edge coverage | `annotate_gfa.py` |
| Estimate copy number | `estimate_copy_numbers.py` |
| Isolate a manageable tangle | `select_tangle.py` |
| Binary optimization | **HUBO**, not QUBO |
| Validate encoding | `validate_hubo.py` |

`annotate_gfa.py` follows the logic of qpg's `tag_gfa_ga.pl`: internal nodes on
an alignment receive full segment coverage, the remaining aligned query length
is distributed over the first/last node, and graph-edge traversal counts are
recorded.

`estimate_copy_numbers.py` follows the logic of qpg's
`tag_gfa_copy_numbers.pl`: estimate a one-copy sequence depth while excluding
obvious cycles/repeats where possible, then fit the depth unit so observed
depths are close to integer multiples.

The original paper also used `pathfinder` for subgraph connection analysis and
depth normalization. The public qpg workflow invokes a pathfinder build with
project-specific options. To keep this starter reproducible, the default
subgraph/copy-number stage is implemented in Python. If you have the same
qpg-compatible `pathfinder` binary, that stage can be swapped in directly while
leaving the Minigraph-Cactus, Giraffe and HUBO stages unchanged.

## HUBO implemented

For a selected graph with `N` biological nodes:

```text
n = ceil(log2(2N))
```

At time `t`, the `n` bits encode the oriented graph state `X_t`.

State IDs are:

```text
2*l     = node l, positive orientation
2*l + 1 = node l, negative orientation
```

The indicator polynomial is

```text
I(X_t=i) = product_k (1 - b_k - x_tk + 2*b_k*x_tk)
```

The generated objective is

```text
H = H_edge + H_copy
```

with

```text
H_edge =
  lambda_edge * sum_t [
      1 - sum_(i,j in E) I(X_t=i) I(X_(t+1)=j)
  ]
```

and

```text
H_copy =
  sum_l [
      sum_t (I(X_t=2l) + I(X_t=2l+1)) - w_l
  ]^2
```

The output is an explicit Boolean polynomial. Terms are automatically reduced
with the binary identity `x^2=x`.

## Inputs

Edit `samples.tsv`:

```text
sample   fasta                         role
REF      data/assemblies/ref.fa        train
TRAIN2   data/assemblies/train2.fa     train
TRAIN3   data/assemblies/train3.fa     train
TRAIN4   data/assemblies/train4.fa     train
TEST     data/assemblies/test.fa       test
```

The held-out sample must **not** be used to build the pangenome.

For an initial demonstration, use a small collection of closely related
haploid assemblies/haplotypes and keep the eventual tangle small. The workflow
uses `--noSplit --permissiveContigFilter` because deliberately small inputs can
otherwise be filtered aggressively by Minigraph-Cactus. Revisit these flags for
full chromosome-scale studies.

## Run

### WSL2 / Ubuntu quick start

This repository is designed to run on **WSL2 Ubuntu on an x86-64 Windows machine**. Keep the clone inside the Linux filesystem (for example `~/projects/`) rather than `/mnt/c/...` because bioinformatics workflows create many intermediate files. Minigraph-Cactus is currently x86-64 oriented.

From a fresh WSL2 terminal:

```bash
git clone https://github.com/aaravshandilya/PangenomeHUBOPipeline-ModifiedPaper.git
cd PangenomeHUBOPipeline-ModifiedPaper
chmod +x setup_wsl.sh
./setup_wsl.sh
conda activate pangenome-hubo
```

Then put your FASTA assemblies in `data/assemblies/`, edit `samples.tsv`, and edit the reference/holdout names in `config.yaml`.

Dry-run the DAG first:

```bash
snakemake -n -p --software-deployment-method conda --cores 4
```

Run the pipeline:

```bash
snakemake -p --software-deployment-method conda --cores 8
```

The first real run can take a while because Snakemake will create the pinned Cactus, vg, and Python environments.

## Why vg indexing is a separate rule

Minigraph-Cactus is used strictly as the **graph builder**. The workflow asks it for the full GFA, then runs `vg autoindex --workflow sr-giraffe` in the pinned vg environment. This avoids coupling the mapping indexes to the version of vg bundled inside the Cactus package and follows the current vg GFA-to-Giraffe indexing workflow.

## Main outputs

```text
results/pangenome/cactus/smallpg.full.gfa.gz
results/pangenome/pangenome.full.gfa
results/pangenome/vg_index.giraffe.gbz

results/annotation/annotated.gfa
results/annotation/node_depths.tsv

results/copy_number/copy_numbers.tsv

results/subgraph/selected.gfa
results/subgraph/selected_copy_numbers.tsv
results/subgraph/summary.json

results/hubo/hubo.json
results/hubo/terms.tsv
results/hubo/variable_map.tsv
results/hubo/metadata.json
results/hubo/validation.json
```

### `hubo.json`

Portable representation:

```json
{
  "format": "binary_polynomial_v1",
  "variables": [
    {"id": 0, "name": "x_0_0", "time": 0, "bit": 0}
  ],
  "terms": [
    {"variables": [], "coefficient": 12.3},
    {"variables": [0], "coefficient": -4.0},
    {"variables": [0, 3, 8], "coefficient": 2.0}
  ],
  "metadata": {}
}
```

An empty variable list is the constant term.

### `terms.tsv`

Human-readable version:

```text
degree  coefficient  variables
0       ...          
1       ...          x_0_0
3       ...          x_0_0,x_0_3,x_1_0
```

## Safety guards

HUBO term counts can grow quickly. `config.yaml` contains:

```yaml
hubo:
  max_variables: 40
  max_terms: 2000000
```

Start small. If the selected graph causes the workflow to stop at this guard,
reduce `subgraph.max_nodes` or explicitly set a smaller walk length.

## Important methodological note

Minigraph-Cactus graph construction is a deliberate change from the original
PGSA paper, which used Minigraph. The downstream optimization problem is still
the same mathematical object required by the HUBO: a weighted oriented graph
`G=(V+ union V-, E, w)`. The graph builder is therefore upstream of the HUBO
encoding rather than part of the HUBO mathematics.
