# In-Silico Genetic Variance Pipeline

This repository prepares variant sets, scores them with AlphaGenome, and aggregates variant-level predictions into downstream gene-level analyses.

## Pipeline Overview

The repository is organized into three stages:

1. **Variant preparation** in `1_variants_prep_nextflow/`
2. **AlphaGenome scoring** in `2_score_variants_with_alphagenome/`
3. **Downstream aggregation** in `3_prepare_results/`

The scoring layer now supports two separate workflows:

- A **legacy batch workflow** for scoring many variants with `score_variants()` plus `GeneMaskLFCScorer`
- A **singleton/haplotype workflow** for comparing a singleton’s baseline score against the same singleton on top of a fine-mapped haplotype background

## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The project expects an AlphaGenome API key in `2_score_variants_with_alphagenome/.env`:

```bash
API_KEY_PERSONAL=your_alphagenome_api_key_here
```

Raw RNA scoring outputs are now treated as protocol-aware downstream data:

- `track_key` is the stable exact RNA track identity reconstructed from AlphaGenome metadata
- `protocol_group` is the grouped RNA protocol label, such as `polyA_plus_rna_seq` or `total_rna_seq`

This keeps distinct RNA protocols from being silently collapsed during aggregation.

## Dataset-Centric Layout

The scoring scripts resolve inputs and outputs through `helpers/path_manager.py`.

Set these before running the scoring scripts:

```bash
export DATASET_ID=dataset_name
export SAMPLE_ID=sample_name
export ROOT_DIR=/path/to/work_root   # optional; falls back to $ROOT_DIR, then $PDC_TMP, then repo-local tmp/
```

This creates a dataset-centric tree:

```text
$ROOT_DIR/experiments/$DATASET_ID/$SAMPLE_ID/
  01_inputs/
  02_chunks/
  03_results/
```

Legacy scoring inputs are discovered from `01_inputs/` using these conventions:

- Variant table: `variants.tsv.gz` or `*_variants.tsv.gz`
- Target BED: `targets.bed` or `*_tss*.bed`
- Optional gene whitelist: `gene_list.tsv` or `ensg_list.txt`

## Stage 1: Variant Preparation

`1_variants_prep_nextflow/` prepares the generic legacy inputs used by the legacy scorer.

Example:

```bash
cd 1_variants_prep_nextflow
nextflow run main.nf \
  --genelist /path/to/gene_list.tsv \
  --mane /path/to/MANE.GRCh38.v1.4.summary.txt \
  --runtime_conf conf/runtime_example.yml
```

If `DATASET_ID` and `SAMPLE_ID` are set, the Nextflow outputs land under the dataset-centric `01_inputs/` directory used by the scorer.

## Stage 2: AlphaGenome Scoring

### Legacy Workflow

The legacy entrypoint remains:

`2_score_variants_with_alphagenome/alphagenome_batch_scoring_batch_window.py`

It still uses:

- `model.score_variants(...)`
- `variant_scorers.GeneMaskLFCScorer(requested_output=RNA_SEQ)`
- chunked `chunk_*.tsv.gz` writing
- optional `STITCH=1` final stitching

Example:

```bash
export DATASET_ID=dataset4
export SAMPLE_ID=background_gnomad
export ROOT_DIR=/scratch/alphagenome

./.venv/bin/python 2_score_variants_with_alphagenome/alphagenome_batch_scoring_batch_window.py
```

Optional overrides:

```bash
./.venv/bin/python 2_score_variants_with_alphagenome/alphagenome_batch_scoring_batch_window.py \
  --variants /path/to/custom_variants.tsv.gz \
  --targets-bed /path/to/custom_targets.bed \
  --gene-list /path/to/gene_list.tsv \
  --batch-size 128 \
  --seq-len-gm 1048576 \
  --stitch
```

Parallel legacy runs still use the existing environment variables:

```bash
export JOB_INDEX=0
export JOB_TOTAL=4
export RAW_START=0
export RAW_END=20000
```

Outputs:

- Chunks: `03_results/../02_chunks/chunk_*.tsv.gz`
- Stitched table: `03_results/$SAMPLE_ID.variants.annotated.tsv.gz`

### Singleton/Haplotype Workflow

This workflow is intentionally separate from the legacy scorer.

It has two steps:

1. Prepare manifests from upstream selected singleton/background parquet files
2. Run the singleton baseline + haplotype scoring entrypoint

#### 1. Prepare Singleton/Haplotype Inputs

Use:

`2_score_variants_with_alphagenome/prepare_singleton_haplotype_inputs.py`

This script consumes **upstream-prepared inputs** from the analysis repository. It does **not** recompute biological selection logic in this repository.

Expected upstream inputs:

- `alphagenome_jobs.parquet`
- `selected_singletons_haplotype_background.parquet`

Example:

```bash
export DATASET_ID=eqtl_singletons
export SAMPLE_ID=singleton_haplotype_v1
export ROOT_DIR=/scratch/alphagenome

./.venv/bin/python 2_score_variants_with_alphagenome/prepare_singleton_haplotype_inputs.py \
  --jobs-parquet /path/to/alphagenome_jobs.parquet \
  --background-parquet /path/to/selected_singletons_haplotype_background.parquet
```

Prepared outputs are written under:

```text
01_inputs/singleton_haplotype/
  singleton_haplotype_jobs.parquet
  singleton_haplotype_jobs.tsv.gz
  singleton_haplotype_background_variants.parquet
  singleton_haplotype_background_variants.tsv.gz
```

The jobs manifest includes:

- `selection_id`, `selection_group`, `group_rank`
- `gene_id`, `gene_symbol`
- singleton `chrom`, `pos`, `ref`, `alt`, `variant_id`
- 1 Mb interval coordinates (`interval_chrom`, `interval_start`, `interval_end`, `interval_width`)
- `background_variants_json`
- normalized upstream summary columns such as `raw_score`, `singleton_af`, and background counts

#### 2. Run Singleton/Haplotype Scoring

Use:

`2_score_variants_with_alphagenome/alphagenome_singleton_haplotype_scoring.py`

This entrypoint runs two explicit branches:

- **Baseline**: `score_variant(interval, singleton_variant, ...)`
- **Haplotype**: build a 1 Mb sequence with background variants applied, predict `background_only`, predict `background_plus_singleton`, then reduce the paired RNA outputs into a track-level scalar

Example:

```bash
export DATASET_ID=eqtl_singletons
export SAMPLE_ID=singleton_haplotype_v1
export ROOT_DIR=/scratch/alphagenome

./.venv/bin/python 2_score_variants_with_alphagenome/alphagenome_singleton_haplotype_scoring.py \
  --jobs-manifest /scratch/alphagenome/experiments/eqtl_singletons/singleton_haplotype_v1/01_inputs/singleton_haplotype/singleton_haplotype_jobs.parquet \
  --reference-fasta /path/to/GRCh38.fa \
  --reference-fai /path/to/GRCh38.fa.fai
```

Useful options:

```bash
./.venv/bin/python 2_score_variants_with_alphagenome/alphagenome_singleton_haplotype_scoring.py \
  --mode baseline \
  --selection-id 2 \
  --selection-id 14
```

```bash
./.venv/bin/python 2_score_variants_with_alphagenome/alphagenome_singleton_haplotype_scoring.py \
  --mode haplotype \
  --ontology-curie UBERON:0011907 \
  --ontology-curie CL:0000515 \
  --ontology-curie CL:0000594
```

By default the new scorer targets:

- Human (`HOMO_SAPIENS`)
- RNA outputs only (`RNA_SEQ`)
- Explicit skeletal-muscle ontologies:
  - `UBERON:0011907`
  - `UBERON:0001134`
  - `CL:0000515`
  - `CL:0000594`

Reference FASTA notes:

- The haplotype branch requires `--reference-fasta`
- If `pysam` is available, indexed compressed FASTA is acceptable
- Without `pysam`, the fallback reader expects an uncompressed FASTA with a `.fai` index

#### New Workflow Outputs

Results are separated from the legacy workflow under:

```text
03_results/singleton_haplotype/
```

Files:

- `baseline_score_variant_raw.parquet` / `.tsv.gz`
- `haplotype_predict_sequence_track_outputs.parquet` / `.tsv.gz`
- `selection_track_comparison.parquet` / `.tsv.gz`
- `selection_track_comparison_wide.parquet` / `.tsv.gz`

Meaning:

- `baseline_score_variant_raw.*` contains the baseline `score_variant()` rows for the selected gene and skeletal-muscle RNA tracks
- `haplotype_predict_sequence_track_outputs.*` contains raw-ish per-track `predict_sequence()` outputs for both `background_only` and `background_plus_singleton`, including track metadata, summary statistics, and serialized track values
- `selection_track_comparison.*` contains long-format comparison rows keyed by `selection_id` and stable `track_id`
- `selection_track_comparison_wide.*` contains one row per `selection_id` + `track_id` with `baseline_scalar`, `haplotype_scalar`, and `haplotype_minus_baseline`

### Reducer Assumption for `predict_sequence()`

`predict_sequence()` returns track predictions, not the legacy `GeneMaskLFCScorer` gene-level delta directly.

To keep that assumption explicit, the new scorer uses:

- `background_only` sequence prediction
- `background_plus_singleton` sequence prediction
- `scalar_score = mean(track_background_plus_singleton - track_background_only)`

This reducer is recorded as:

`reducer_name = mean_alt_minus_background_full_track`

The comparison tables also include `scalar_score_center_bin_delta` so downstream analysis can inspect a more local effect summary without changing the primary reducer silently.

## Stage 3: Aggregation & Analysis

The downstream aggregation modules consume scored variant tables and compute metrics such as predicted genetic variance:

$$ V_{G,pred} = \sum_{i \in \text{variants}} 2 \cdot p_i \cdot (1 - p_i) \cdot \beta_i^2 $$

Where:

- \( p_i \) is allele frequency
- \( \beta_i \) is the predicted effect size (`raw_score` in the legacy workflow)

## Notes

- The legacy workflow was kept separate from the singleton/haplotype workflow to avoid silently changing existing runs.
- The singleton/haplotype prep step consumes upstream selected manifests from the analysis repository; this repository only normalizes them for scoring.
- `2_score_variants_with_alphagenome/batch_stability_scorer.py` now also resolves its legacy inputs from the same dataset-centric layout instead of hardcoded absolute paths.
- Downstream aggregation is protocol-aware: exact RNA rows are preserved per `track_key`, while gene-level summaries can still be grouped by `protocol_group`.
- If an older variant parquet already dropped RNA track identity, rebuild it from raw chunk TSVs or rerun scoring before aggregating.
