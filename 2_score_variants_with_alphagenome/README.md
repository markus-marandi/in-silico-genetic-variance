# AlphaGenome Scoring Refactor

This directory now supports two workflow families:

1. Legacy batch scoring for the original gene-set pipeline
2. Singleton-focused scoring for the eQTL / haplotype-background experiment

The refactor keeps the legacy workflow available, removes hardcoded input paths, and adds an explicit prep-plus-score path for singleton baseline and singleton haplotype analyses.

## What Changed From the Old Single-Script Workflow

### Before

The old workflow was effectively one batch script plus one Python entrypoint:

- `alphagenome_batch_scoring_batch_window.py`
- hardcoded `VAR_TSV`, `TSS_BED`, and `GENE_LIST_PATH`
- optional hardcoded defaults for `DATASET_ID` / `SAMPLE_ID`
- one scoring style: reference-based `score_variants(...)`
- chunked execution driven by `JOB_INDEX`, `JOB_TOTAL`, `RAW_START`, and `RAW_END`

Typical submission style looked like:

```bash
sbatch --export=ENV_FILE=${ENV_PATH},API_KEY_VAR=FIRST_KEY,JOB_INDEX=0,JOB_TOTAL=3,DATASET_ID=dataset4,SAMPLE_ID=background run_alphagenome_batch_window.sh
```

### After

The refactor keeps that legacy path, but changes how inputs are resolved and adds a separate singleton experiment path:

- `alphagenome_batch_scoring_batch_window.py`
  - still handles legacy batch scoring
  - now resolves inputs from the dataset-centric layout or explicit CLI overrides
  - no Python edits needed for per-run paths
- `prepare_singleton_haplotype_inputs.py`
  - converts upstream singleton/background parquet inputs into scorer-ready manifests
- `alphagenome_singleton_haplotype_scoring.py`
  - `--mode baseline` runs singleton-against-reference scoring with `score_variant(...)`
  - `--mode haplotype` builds a modified 1 Mb sequence and runs `predict_sequence(...)`
  - `--mode both` runs both branches in one pass

## Preferred Architecture

The intended architecture is:

- legacy batch scoring stays available and stable
- a prep step creates singleton / haplotype manifests under `01_inputs/singleton_haplotype/`
- a separate singleton scorer handles the new experiment explicitly
- baseline and haplotype branches stay conceptually separate even though the current singleton scorer exposes them through one `--mode`

Today, the repo implements this as separate Python entrypoints plus a prep step. If you want one future-state SLURM launcher, the recommended wrapper convention is:

- `WORKFLOW=legacy_batch`
- `WORKFLOW=singleton_baseline`
- `WORKFLOW=singleton_haplotype`

with the wrapper dispatching to the correct Python command. The Python layer itself does not require `WORKFLOW`; it uses separate scripts plus `--mode`.

## Dataset-Centric Layout

All workflows use `helpers/path_manager.py` and the shared `ProjectLayout` conventions:

```text
$ROOT_DIR/experiments/$DATASET_ID/$SAMPLE_ID/
  01_inputs/
  02_chunks/
  03_results/
```

Legacy scoring resolves its default inputs from `01_inputs/`:

- `variants.tsv.gz` or another matching `*_variants.tsv.gz`
- `targets.bed` or another matching `*_tss*.bed`
- optional `gene_list.tsv` or `ensg_list.txt`

Singleton prep writes normalized manifests to:

```text
01_inputs/singleton_haplotype/
  singleton_haplotype_jobs.parquet
  singleton_haplotype_jobs.tsv.gz
  singleton_haplotype_background_variants.parquet
  singleton_haplotype_background_variants.tsv.gz
```

Legacy outputs remain split between:

- `02_chunks/chunk_*.tsv.gz`
- `03_results/${SAMPLE_ID}.variants.annotated.tsv.gz` when stitching is enabled

Singleton outputs are isolated under:

```text
03_results/singleton_haplotype/
```

including:

- `baseline_score_variant_raw.parquet` / `.tsv.gz`
- `haplotype_predict_sequence_track_outputs.parquet` / `.tsv.gz`
- `selection_track_comparison.parquet` / `.tsv.gz`
- `selection_track_comparison_wide.parquet` / `.tsv.gz`

## Legacy Usage Still Works

The legacy scorer remains the same operational entrypoint:

- same script name: `alphagenome_batch_scoring_batch_window.py`
- same chunking model
- same `score_variants(...)` plus `GeneMaskLFCScorer`
- same `JOB_INDEX` / `JOB_TOTAL` fan-out pattern
- same optional `STITCH=1` behavior

What changed is input resolution:

- old behavior: edit hardcoded paths in Python
- new behavior: place files under the dataset-centric `01_inputs/` directory or pass explicit `--variants`, `--targets-bed`, and `--gene-list` overrides

This means old batch usage stays conceptually the same, but path selection moves out of the script body and into the run environment.

## Current SLURM UX Baseline

The baseline launcher style is the pair of templates:

- [`run_alphagenome_batch_window_example.sh`](./run_alphagenome_batch_window_example.sh) for legacy batch scoring
- [`run_alphagenome_singleton_haplotype_example.sh`](./run_alphagenome_singleton_haplotype_example.sh) for singleton baseline / haplotype scoring

They do four important things:

- loads your modules and Python environment
- loads a `.env` file via `ENV_FILE`
- selects the API key by name via `API_KEY_VAR` or directly via `API_KEY_OVERRIDE`
- passes `DATASET_ID`, `SAMPLE_ID`, `JOB_INDEX`, and `JOB_TOTAL` through `sbatch --export=...`

The template is intentionally incomplete for site-specific details such as `ENV_DIR`, `cd`, and the final Python command. The examples below keep that exact submission style and only change the target command or add a small number of new exports.

## Environment Variables

### Common Runtime Variables

These are the main variables used by the current launcher pattern and Python entrypoints:

| Variable | Required | Used by | Meaning |
| --- | --- | --- | --- |
| `DATASET_ID` | yes | all workflows unless passed by CLI | logical dataset group |
| `SAMPLE_ID` | yes | all workflows unless passed by CLI | run / sample name |
| `ROOT_DIR` | no | all workflows | base experiment root |
| `PDC_TMP` | no | all workflows | fallback root when `ROOT_DIR` is unset |
| `ENV_FILE` | launcher-level | SLURM template | `.env` file to source before choosing the key |
| `API_KEY_VAR` | launcher-level | SLURM template | env var name to read from `.env` |
| `API_KEY_OVERRIDE` | launcher-level | SLURM template | direct API key override |
| `API_KEY_PERSONAL` | yes after env loading | Python scripts | effective AlphaGenome API key env var |

### Legacy Batch Variables

| Variable | Required | Used by | Meaning |
| --- | --- | --- | --- |
| `JOB_INDEX` | no | legacy batch | zero-based job index, default `0` |
| `JOB_TOTAL` | no | legacy batch | total job fan-out, default `1` |
| `RAW_START` | no | legacy batch | manual start index override before batch alignment |
| `RAW_END` | no | legacy batch | manual end index override before batch alignment |
| `BATCH` | no | legacy batch | fallback default for `--batch-size` |
| `SEQ_LEN_GM` | no | legacy batch | fallback default for `--seq-len-gm` |
| `STITCH` | no | legacy batch | set `1` to stitch chunk outputs |

### Launcher Template Variables

These come from the current shell template rather than the Python code:

| Variable | Required | Meaning |
| --- | --- | --- |
| `ENV_DIR` | yes in the template | conda or venv root used to find Python |
| `PYTHON` | derived | `${ENV_DIR}/bin/python` in the template |

### Recommended Wrapper Variables for a Single Future-State Launcher

These are recommended naming conventions if you keep one SLURM wrapper for all workflows:

| Variable | Required | Meaning |
| --- | --- | --- |
| `WORKFLOW` | recommended | `legacy_batch`, `singleton_baseline`, or `singleton_haplotype` |
| `REFERENCE_FASTA` | haplotype only | reference FASTA passed through to `--reference-fasta` |
| `REFERENCE_FAI` | optional | FASTA index passed through to `--reference-fai` |

`WORKFLOW`, `REFERENCE_FASTA`, and `REFERENCE_FAI` are wrapper-level conventions. The current Python scripts use CLI arguments rather than reading those variables directly.

## CLI Arguments

### `alphagenome_batch_scoring_batch_window.py`

Legacy batch scorer:

- `--dataset-id`
- `--sample-id`
- `--root-dir`
- `--variants`
- `--targets-bed`
- `--gene-list`
- `--env-file`
- `--batch-size`
- `--seq-len-gm`
- `--stitch`

### `prepare_singleton_haplotype_inputs.py`

Singleton manifest prep:

- `--dataset-id`
- `--sample-id`
- `--root-dir`
- `--jobs-parquet` required
- `--background-parquet` required
- `--seq-len`
- `--output-dir`

### `alphagenome_singleton_haplotype_scoring.py`

Singleton scoring:

- `--dataset-id`
- `--sample-id`
- `--root-dir`
- `--jobs-manifest`
- `--background-manifest`
- `--reference-fasta` required for `--mode haplotype` and `--mode both`
- `--reference-fai`
- `--mode {both,baseline,haplotype}`
- `--env-file`
- `--ontology-curie` repeatable
- `--selection-id` repeatable
- `--max-selections`
- `--output-dir`

## How To Run Each Workflow

### 1. Legacy Batch Scoring

This stays closest to the current submission style.

If you keep your existing launcher name:

```bash
sbatch --export=ENV_FILE=${ENV_PATH},API_KEY_VAR=API_KEY_PERSONAL_1,JOB_INDEX=0,JOB_TOTAL=3,DATASET_ID=dataset4,SAMPLE_ID=background run_alphagenome_batch_window.sh
```

The launcher should call:

```bash
PYTHONUNBUFFERED=1 "${PYTHON}" \
  2_score_variants_with_alphagenome/alphagenome_batch_scoring_batch_window.py
```

Legacy mode still supports explicit path overrides when needed:

```bash
PYTHONUNBUFFERED=1 "${PYTHON}" \
  2_score_variants_with_alphagenome/alphagenome_batch_scoring_batch_window.py \
  --variants /path/to/custom_variants.tsv.gz \
  --targets-bed /path/to/custom_targets.bed \
  --gene-list /path/to/gene_list.tsv \
  --batch-size 128 \
  --seq-len-gm 1048576 \
  --stitch
```

### 2. Prepare Singleton / Haplotype Manifests

This is the new prep step. It consumes upstream-prepared parquets and writes scorer-ready manifests under `01_inputs/singleton_haplotype/`.

Direct invocation:

```bash
PYTHONUNBUFFERED=1 "${PYTHON}" \
  2_score_variants_with_alphagenome/prepare_singleton_haplotype_inputs.py \
  --dataset-id eqtl_singletons \
  --sample-id singleton_haplotype_v1 \
  --root-dir /scratch/alphagenome \
  --jobs-parquet /path/to/alphagenome_jobs.parquet \
  --background-parquet /path/to/selected_singletons_haplotype_background.parquet
```

This prep step should be run before any singleton baseline or haplotype scoring.

### 3. Singleton Baseline Scoring

This branch scores each selected singleton against reference using `score_variant(...)`.

Current Python command:

```bash
PYTHONUNBUFFERED=1 "${PYTHON}" \
  2_score_variants_with_alphagenome/alphagenome_singleton_haplotype_scoring.py \
  --dataset-id eqtl_singletons \
  --sample-id singleton_haplotype_v1 \
  --root-dir /scratch/alphagenome \
  --mode baseline
```

If you want this to feel like the current `sbatch --export=...` UX, the recommended wrapper convention is:

```bash
sbatch --export=ENV_FILE=${ENV_PATH},API_KEY_VAR=API_KEY_PERSONAL_1,WORKFLOW=singleton_baseline,DATASET_ID=eqtl_singletons,SAMPLE_ID=singleton_haplotype_v1 run_alphagenome_batch_window.sh
```

There is also a dedicated singleton launcher template:

```bash
sbatch --export=ENV_FILE=${ENV_PATH},API_KEY_VAR=API_KEY_PERSONAL_1,MODE=baseline,DATASET_ID=eqtl_singletons,SAMPLE_ID=singleton_haplotype_v1 run_alphagenome_singleton_haplotype_example.sh
```

Difference from legacy:

- same submission shape
- no hardcoded input paths
- uses prepared singleton manifests instead of generic variant tables
- usually no `JOB_INDEX` / `JOB_TOTAL` fan-out is needed
- current implementation does not consume `JOB_INDEX` / `JOB_TOTAL` in singleton mode

### 4. Singleton Haplotype Scoring

This branch builds a modified background sequence and runs `predict_sequence(...)` on:

- `background_only`
- `background_plus_singleton`

Current Python command:

```bash
PYTHONUNBUFFERED=1 "${PYTHON}" \
  2_score_variants_with_alphagenome/alphagenome_singleton_haplotype_scoring.py \
  --dataset-id eqtl_singletons \
  --sample-id singleton_haplotype_v1 \
  --root-dir /scratch/alphagenome \
  --mode haplotype \
  --reference-fasta /path/to/GRCh38.fa \
  --reference-fai /path/to/GRCh38.fa.fai
```

Recommended single-wrapper submission style:

```bash
sbatch --export=ENV_FILE=${ENV_PATH},API_KEY_VAR=API_KEY_PERSONAL_1,WORKFLOW=singleton_haplotype,DATASET_ID=eqtl_singletons,SAMPLE_ID=singleton_haplotype_v1,REFERENCE_FASTA=/path/to/GRCh38.fa,REFERENCE_FAI=/path/to/GRCh38.fa.fai run_alphagenome_batch_window.sh
```

Dedicated singleton launcher style:

```bash
sbatch --export=ENV_FILE=${ENV_PATH},API_KEY_VAR=API_KEY_PERSONAL_1,MODE=haplotype,DATASET_ID=eqtl_singletons,SAMPLE_ID=singleton_haplotype_v1,REFERENCE_FASTA=/path/to/GRCh38.fa,REFERENCE_FAI=/path/to/GRCh38.fa.fai run_alphagenome_singleton_haplotype_example.sh
```

Difference from legacy:

- same `sbatch --export=...` feel
- requires the prep step first
- requires a reference FASTA for sequence construction
- uses `predict_sequence(...)` instead of the legacy batch `score_variants(...)`
- produces track-level comparison outputs instead of chunked variant tables

## Recommended Single-Launcher Dispatch

If you want to preserve one familiar shell wrapper, keep the current environment-loading logic and dispatch the final Python call by workflow:

```bash
case "${WORKFLOW:-legacy_batch}" in
  legacy_batch)
    PYTHONUNBUFFERED=1 "${PYTHON}" \
      2_score_variants_with_alphagenome/alphagenome_batch_scoring_batch_window.py
    ;;
  singleton_baseline)
    PYTHONUNBUFFERED=1 "${PYTHON}" \
      2_score_variants_with_alphagenome/alphagenome_singleton_haplotype_scoring.py \
      --mode baseline
    ;;
  singleton_haplotype)
    if [ -n "${REFERENCE_FAI:-}" ]; then
      PYTHONUNBUFFERED=1 "${PYTHON}" \
        2_score_variants_with_alphagenome/alphagenome_singleton_haplotype_scoring.py \
        --mode haplotype \
        --reference-fasta "${REFERENCE_FASTA}" \
        --reference-fai "${REFERENCE_FAI}"
    else
      PYTHONUNBUFFERED=1 "${PYTHON}" \
        2_score_variants_with_alphagenome/alphagenome_singleton_haplotype_scoring.py \
        --mode haplotype \
        --reference-fasta "${REFERENCE_FASTA}"
    fi
    ;;
esac
```

That keeps the user-facing SLURM UX close to the current style while making the new branches explicit.

## Output Semantics

### Legacy Batch

Legacy outputs remain chunked variant-level tables, with optional final stitching into:

- `03_results/${SAMPLE_ID}.variants.annotated.tsv.gz`

### Singleton Baseline

Baseline outputs are raw `score_variant()` rows for the selected singleton and target gene, filtered to skeletal-muscle RNA tracks.

Important fields include:

- `selection_id`
- `selection_group`
- `gene_id`
- `gene_symbol`
- `singleton_variant_id`
- `track_id` as a protocol-aware stable track identifier
- `raw_score`
- `scalar_score`
- `reducer_name = GeneMaskLFCScorer.raw_score`

### Singleton Haplotype

Haplotype outputs contain raw-ish per-track `predict_sequence()` summaries for both:

- `background_only`
- `background_plus_singleton`

and reduced comparison tables with:

- `baseline_scalar`
- `haplotype_scalar`
- `haplotype_minus_baseline`
- `scalar_score_center_bin_delta`

The current haplotype reducer is explicitly documented as:

- `reducer_name = mean_alt_minus_background_full_track`

## Migration Notes From the Old Hardcoded-Path Workflow

1. Stop editing Python source to change `VAR_TSV`, `TSS_BED`, or `GENE_LIST_PATH`.
2. Put legacy inputs under the dataset-centric `01_inputs/` directory, or pass explicit CLI overrides.
3. Keep using the same `sbatch --export=...` pattern for legacy batch runs.
4. Add a prep step before any singleton scoring. The singleton scorer expects prepared manifests, not generic whole-gene variant TSVs.
5. Pass `--reference-fasta` for haplotype mode. Baseline mode does not need it.
6. Do not expect singleton scoring to reuse legacy anchor-row logic. The singleton workflow consumes prepared intervals directly.
7. Treat `WORKFLOW` as a wrapper convenience, not a Python requirement.
8. Expect singleton outputs under `03_results/singleton_haplotype/`, not alongside legacy chunk files.

## Practical Summary

- Use `alphagenome_batch_scoring_batch_window.py` for the legacy gene-set batch workflow.
- Use `prepare_singleton_haplotype_inputs.py` to materialize singleton manifests.
- Use `alphagenome_singleton_haplotype_scoring.py --mode baseline` for singleton-against-reference scoring.
- Use `alphagenome_singleton_haplotype_scoring.py --mode haplotype` for background-sequence plus singleton scoring.
- If you want one wrapper script, keep your current `sbatch --export=...` UX and add a thin `WORKFLOW` dispatch layer.
