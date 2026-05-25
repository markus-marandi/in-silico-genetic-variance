# In-Silico Genetic Variance Pipeline

This repository contains the production pipeline used to prepare population and synthetic variant sets, score them with AlphaGenome, and aggregate variant-level expression effect predictions into gene-level estimates of predicted cis-regulatory genetic variance (`V_G`).

The pipeline supports the thesis **Deep-Learning Prediction of Human Gene Expression Variability from Genetic Variation** (Markus Marandi, KTH Royal Institute of Technology, 2026). Its companion analysis repository is [`in-silico-vg-analysis`](https://github.com/markus-marandi/in-silico-vg-analysis).

![Computational workflow for variant generation and scoring](<docs/figs/Figure 1.9.1: Computational workflow for variant generation and scoring.png>)

## Scientific Aim

Many haploinsufficient genes are sensitive to reduced gene dosage, but existing gene-constraint metrics focus mainly on protein-coding intolerance, conservation, structural variation, or static epigenomic annotations. This project asks whether sequence-based deep learning can estimate a complementary gene-level measure: the depletion of predicted cis-regulatory expression variance.

The central quantity is predicted genetic variance:

```text
V_G = sum_i 2 p_i (1 - p_i) beta_i^2
```

where `p_i` is the allele frequency of variant `i`, and `beta_i` is the AlphaGenome-predicted RNA-seq effect size for that variant. In public-facing documentation, this quantity should be described as an effect size, `beta`, or `Delta`; `raw_score` is retained only as an internal column name.

## Thesis Workflow

The thesis workflow proceeds in three stages:

1. **Variant preparation**: define MANE Select v1.4 transcription start site windows on GRCh38, extract observed gnomAD v4.1 variants, and generate context-matched synthetic SNVs.
2. **AlphaGenome scoring**: score observed and synthetic variants using RNA-seq outputs and exon-mask aggregation to obtain per-variant expression effect sizes.
3. **Aggregation**: deduplicate protocol-aware variant outputs, assign allele frequencies, compute gene-level `V_G`, and export variant- and gene-level tables for downstream analysis.

The primary thesis analysis used 316 ClinGen haploinsufficient genes and 349 background coding genes after transcript, expression, and promoter-overlap filtering. Variants were restricted to +/-10 kb around each gene's TSS, while AlphaGenome used its 1 Mb sequence context during inference.

## Repository Structure

```text
in-silico-genetic-variance/
├── 1_variants_prep_nextflow/        # gnomAD extraction and synthetic SNV generation
├── 2_score_variants_with_alphagenome/ # AlphaGenome scoring entrypoints
├── 3_prepare_results/               # variant normalization and gene-level aggregation
├── docs/                            # research notes, public summaries, and workflow figures
├── helpers/                         # shared path-management helpers
├── extras/                          # older or auxiliary utilities
└── requirements.txt
```

## Pipeline Stages

### 1. Variant Preparation

`1_variants_prep_nextflow/` builds TSS-centered target intervals, extracts observed variants from gnomAD v4.1, and can generate synthetic SNVs absent from gnomAD. The synthetic set is constructed by enumerating possible SNVs within each target window, excluding observed population variants by anti-join, matching mutation class and CpG methylation context, and downsampling per gene to match observed variant counts.

### 2. AlphaGenome Scoring

`2_score_variants_with_alphagenome/` contains the scoring layer. The thesis batch workflow scores many observed or synthetic variants with AlphaGenome RNA-seq outputs and gene-level exon-mask aggregation. A newer singleton/haplotype workflow is kept separate because it answers a different, exploratory question about population-background effects around selected singleton variants.

### 3. Downstream Aggregation

`3_prepare_results/` converts scored chunks into analysis-ready variant and gene tables. It preserves exact RNA track identity through `track_key`, groups related RNA protocols through `protocol_group`, and computes gene-level quantities such as total `V_G`, allele-frequency partitions, and regional TSS-window contributions.

For biological observed-versus-null comparisons:

- observed gnomAD datasets use `vg_predicted`, `vg_common`, `vg_rare`, and regional `vg_*` columns;
- synthetic/null datasets use `vg_predicted_perm`, `vg_common_perm`, `vg_rare_perm`, and regional `vg_*_perm` columns;
- permutation/sanity columns in real files are quality-control artifacts, not the biological null model.

## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

AlphaGenome scoring expects an API key in `2_score_variants_with_alphagenome/.env`:

```bash
API_KEY_PERSONAL=your_alphagenome_api_key_here
```

Stage 1 also requires a Hail-compatible environment for gnomAD extraction. See `1_variants_prep_nextflow/README.md` for Nextflow and runtime configuration details.

## Usage Sketch

Prepare observed variants:

```bash
cd 1_variants_prep_nextflow
nextflow run main.nf \
  --genelist /path/to/gene_list.tsv \
  --mane /path/to/MANE.GRCh38.v1.4.summary.txt \
  --runtime_conf conf/runtime_example.yml
```

Score variants with the batch AlphaGenome workflow:

```bash
export DATASET_ID=dataset4
export SAMPLE_ID=background_gnomad
export ROOT_DIR=/scratch/alphagenome

python 2_score_variants_with_alphagenome/alphagenome_batch_scoring_batch_window.py
```

Aggregate scored variants:

```bash
python 3_prepare_results/pipeline_runner.py \
  --variants-parquet /path/to/scored_variants.parquet \
  --variant-out /path/to/variants_dedup.parquet \
  --gene-out /path/to/gene_level_metrics.parquet \
  --deduplicate \
  --permute-af \
  --calc-ci
```

The exact production commands depend on the dataset, HPC filesystem, and whether the input is observed gnomAD or synthetic/null data.

## Interpretation Boundaries

Predicted `V_G` is a gene-level research screening metric. It is not a standalone clinical score, and it should not be used as independent evidence of pathogenicity without functional or population-genetic validation.

The thesis found statistically significant but modest gene-level concordance between predicted `V_G` and GTEx eQTL-derived `V_G`; per-variant eQTL concordance was not supported in the matched set. Synthetic-null comparisons should also be interpreted cautiously because local sequence composition, especially CpG-rich promoter architecture, can influence both AlphaGenome effect sizes and null expectations.

## Acknowledgements and Funding

This work was supervised by Philipp Rentzsch and Tuuli Lappalainen at the Lappalainen Lab, KTH / Science for Life Laboratory.

Computational resources were provided by the National Academic Infrastructure for Science and Information Systems (NAISS) at PDC Center for High Performance Computing, KTH Royal Institute of Technology, allocation `NAISS 2024/6-322`, and at the National Supercomputer Centre (NSC), Linköping University, allocation `Berzelius-2025-176`.

This work was supported by a Wallenberg Scholar award 2024 to T. Lappalainen from the Knut and Alice Wallenberg Foundation, grant `KAW 2023.0337`, "Functional architecture of genetic disease risk using natural variation and experimental perturbations", and by a Göran Gustavsson award 2023 to T. Lappalainen.

## Citation

If you use this repository, please cite the repository and the associated thesis. Repository metadata is provided in `CITATION.cff`.
