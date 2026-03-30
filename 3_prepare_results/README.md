## Single Entrypoint

The main driver here is `pipeline_runner.py`. It handles two primary modes:

1. **Raw Processing:** Stitching raw chunk files  Variants  Genes.
2. **Post-Processing:** Loading existing Variant Parquets  Deduplication/Permutation  Genes.

### Environment Variables

* `DATASET_ID`: Logical dataset grouping (e.g., `dataset5_null`).
* `SAMPLE_ID`: Run-specific sample name (e.g., `background_ISM`).
* `ROOT_DIR` (Optional): Base scratch root (defaults to `$PDC_TMP` or `/cfs/...`).

---

## Key Features & Modules

* **`stitcher.py`**: Lazy scans `chunk_*.tsv.gz`, normalizes IDs, and applies optional gene whitelisting.
* **`normalizer.py`**: Standardizes columns (`CHROM`, `POS`, `REF`, `ALT`) and ensures method-friendly tagging.
* **`annotator.py`**: Joins Allele Frequency (AF) from inputs or gnomAD.
* **`deduplicator.py` (New)**:
* **Protocol-Aware Deduplication**: Removes duplicate variants only within the same `gene_id`, `variant_id`, and exact `track_key`.
* **Deterministic Tie-Breaking**: Sorts by `abs_score` (desc) then `variant_id` (asc) to ensure 100% reproducibility.


* **`downsampler.py`**:
* **Rigid Downsampling**: For synthetic (ISM/Null) datasets, randomly downsamples variants to match the exact count of variants found in a "Real" reference dataset, gene-by-gene. Uses a fixed seed (`42`).


* **`permuted_af.py`**:
* **AF Permutation**: "Borrows" the AF distribution from a real dataset and randomly assigns it to synthetic variants (preserving the per-gene AF distribution).


* **`aggregator.py`**:
* **Gene Counts & Stats**: Aggregates variance (), spatial bins (Promoter/Upstream), and enrichment scores.
* **Protocol-Aware Long Format**: Gene-level outputs keep separate rows per `gene_id` and `protocol_group`.
* **Dual Metrics**: Calculates `vg_predicted` (using real AF) and `vg_predicted_perm` (using permuted AF) simultaneously.
* **Confidence Intervals**: [Optional] Runs a Monte Carlo simulation (1,000 iterations) during aggregation to calculate the 5th and 95th percentiles for .



---

## Usage

### 1. Standard Run (Chunks  Output)

Process raw scoring chunks into final outputs.

```bash
python pipeline_runner.py \
  --chunks-dir /cfs/.../02_chunks \
  --variants-af /cfs/.../01_inputs/variants.tsv.gz \
  --gnomad-af /cfs/.../gnomad.parquet \
  --gene-list /path/to/gene_list.tsv

```

### 2. Post-Processing / Migration (Existing Parquet  Output)

Use this mode to apply deduplication, downsampling, or AF permutation to existing files.

#### Flags Reference

| Flag | Description |
| --- | --- |
| `--variants-parquet` | Input variant parquet file to process. |
| `--deduplicate` | Triggers rigid deduplication (and downsampling if ISM). |
| `--real-reference` | Path to "Real" dataset. Required for ISM downsampling or AF permutation. |
| `--permute-af` | Generates `perm_AF` column by sampling from the reference pool. |
| `--calc-ci` | Runs Monte Carlo simulation to compute  confidence intervals (p5/p95). |
| `--gene-list` | Optional whitelist to filter variants *before* processing. |

---

## Workflows: Real vs. Synthetic

These 4 commands represent the standard workflow to generate the **Background** and **ClinGen** datasets for "Observed vs Simulated" analysis.

### 1. Background Real (GnomAD)

*Goal: Clean, deduplicate, and create a sanity-check permutation (shuffling AFs within genes).*

```bash
python pipeline_runner.py \
  --variants-parquet /cfs/.../background_variants_20260102.parquet \
  --variant-out /cfs/.../Background_Gnomad_variants_dedup_perm.parquet \
  --gene-out /cfs/.../Background_Gnomad_genes.parquet \
  --deduplicate \
  --permute-af \
  --calc-ci

```

### 2. ClinGen Real (GnomAD)

*Goal: Filter to specific ClinGen genes, deduplicate, and create sanity-check permutation.*

```bash
python pipeline_runner.py \
  --variants-parquet /cfs/.../clingen_alphagenome_scores.backfilled.parquet \
  --gene-list /cfs/.../ClinGen_gene_curation_list_GRCh38.ensg.txt \
  --variant-out /cfs/.../ClinGen_HI_Gnomad_variants_dedup.parquet \
  --gene-out /cfs/.../ClinGen_HI_Gnomad_genes.parquet \
  --deduplicate \
  --permute-af \
  --calc-ci

```

### 3. Background Synthetic (Null/ISM)

*Goal: Downsample synthetic variants to match Real Background counts, then borrow Real AFs.*

```bash
python pipeline_runner.py \
  --variants-parquet /cfs/.../dataset5_Background_NULL_variant_level_summary.parquet \
  --real-reference /cfs/.../Background_Gnomad_variants_dedup_perm.parquet \
  --variant-out /cfs/.../Background_Synth_variants_downsampled_perm.parquet \
  --gene-out /cfs/.../Background_Synth_genes.parquet \
  --deduplicate \
  --permute-af \
  --calc-ci

```

### 4. ClinGen Synthetic (Null/ISM)

*Goal: Downsample synthetic variants to match Real ClinGen counts, then borrow Real AFs.*

```bash
python pipeline_runner.py \
  --variants-parquet /cfs/.../dataset5_ClinGen_NULL_variant_level_summary.parquet \
  --real-reference /cfs/.../ClinGen_HI_Gnomad_variants_dedup.parquet \
  --variant-out /cfs/.../ClinGen_HI_Synth_variants_downsampled_perm.parquet \
  --gene-out /cfs/.../ClinGen_HI_Synth_genes.parquet \
  --deduplicate \
  --permute-af \
  --calc-ci

```

---

## Output Columns (Gene Level)

The aggregator produces comprehensive gene-level metrics. Key columns include:

* **`protocol_group`**: protocol-level RNA grouping with values `polyA_plus_rna_seq`, `total_rna_seq`, or `other`.
* **`n_track_keys`**: number of exact AlphaGenome RNA tracks contributing to that gene/protocol row.
* **`vg_predicted`**: Genetic Variance calculated using the standard `AF` column.
* **`vg_predicted_perm`**: Genetic Variance calculated using the `perm_AF` column (if `--permute-af` was used).
* **`vg_perm_mean`, `vg_perm_p05`, `vg_perm_p95**`: Monte Carlo statistics for error bars (if `--calc-ci` was used).
* **`n_variants`**: Final count of variants (after deduplication/downsampling).
* **`mean_abs_effect`**: Mean absolute score of variants in the gene.
* **`mean_abs_promoter` / `mean_abs_gene_body**`: Spatial scoring breakdown.

## Protocol-Aware Variant Identity

Variant-level RNA outputs now keep two explicit identity columns:

* **`protocol_group`**: clean grouping for protocol-level summaries.
* **`track_key`**: stable exact track identifier built from AlphaGenome track metadata such as `Assay title` / `Assay_title`, `track_name`, `ontology_curie`, `biosample_name`, and `gtex_tissue`.

This means the same `gene_id` + `variant_id` can legitimately appear multiple times in processed outputs when distinct RNA tracks or protocols are present. `polyA_plus_rna_seq` and `total_rna_seq` are no longer collapsed together.

## Rerun Boundary

* rerun only `3_prepare_results` if your variant parquet or rebuilt chunk-to-parquet step still preserves the raw track metadata needed to derive `track_key`
* rerun scoring if the existing variant parquet already dropped track metadata and the prepare-results stage can no longer reconstruct exact track identity


----------


python pipeline_runner.py \
  --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/background_variants_20260102.parquet \
  --variant-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/Background_Gnomad_variants_dedup_perm.parquet \
  --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/Background_Gnomad_genes.parquet \
  --deduplicate \
  --permute-af


python -c "import polars as pl; df = pl.read_parquet('/cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/Background_Gnomad_variants_dedup_perm.parquet'); print(f'Total Rows: {len(df)}'); print(f'Unique Genes: {df[\"gene_id\"].n_unique()}')"

python -c "import polars as pl; df = pl.read_parquet('/cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/Background_Gnomad_genes.parquet'); print(f'Total Rows: {len(df)}'); print(f'Unique Genes: {df[\"gene_id\"].n_unique()}')"

Total Rows: 1999142
Unique Genes: 349

2. ClinGen Real (Filter + Cleanup)

python pipeline_runner.py \
  --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/clingen_alphagenome_scores_all_aggs_variantids_long.backfilled.parquet \
  --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/01_inputs/ClinGen_gene_curation_list_GRCh38.ensg.txt \
  --variant-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/ClinGen_HI_Gnomad_variants_dedup_26012026.parquet \
  --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/ClinGen_HI_Gnomad_genes_26012026.parquet \
  --deduplicate

python -c "import polars as pl; df = pl.read_parquet('/cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/ClinGen_HI_Gnomad_variants_dedup_26012026.parquet'); print(f'Total Rows: {len(df)}'); print(f'Unique Genes: {df[\"gene_id\"].n_unique()}')"
Total Rows: 1743183
Unique Genes: 316
python -c "import polars as pl; df = pl.read_parquet('/cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/ClinGen_HI_Gnomad_genes_26012026.parquet'); print(f'Total Rows: {len(df)}'); print(f'Unique Genes: {df[\"gene_id\"].n_unique()}')"


3. Background Synth (Downsample + Permute)


python pipeline_runner.py \
  --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/03_results/dataset5_Background_NULL_variant_level_summary.parquet \
  --real-reference /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/Background_Gnomad_variants_dedup_perm_26012026.parquet \
  --variant-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/03_results/Background_Synth_variants_downsampled_perm_26012026.parquet \
  --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/03_results/Background_Synth_genes.parquet \
  --deduplicate \
  --permute-af

python -c "import polars as pl; df = pl.read_parquet('/cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/03_results/Background_Synth_variants_downsampled_perm_26012026.parquet'); print(f'Total Rows: {len(df)}'); print(f'Unique Genes: {df[\"gene_id\"].n_unique()}')"

python -c "import polars as pl; df = pl.read_parquet('/cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/03_results/Background_Synth_genes.parquet'); print(f'Total Rows: {len(df)}'); print(f'Unique Genes: {df[\"gene_id\"].n_unique()}')"

  downsample output: (1987481, 34)
349 genes


4. ClinGen Synth (Downsample + Permute)

python pipeline_runner.py \
  --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/03_results/dataset5_ClinGen_NULL_variant_level_summary.parquet \
  --real-reference /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/ClinGen_HI_Gnomad_variants_dedup_26012026.parquet \
  --variant-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/03_results/ClinGen_HI_Synth_variants_downsampled_perm_26012026.parquet \
  --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/03_results/ClinGen_HI_Synth_genes_26012026.parquet \
  --deduplicate \
  --permute-af
