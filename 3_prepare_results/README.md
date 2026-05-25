# Prepare Results

This stage converts AlphaGenome scoring outputs into analysis-ready variant- and gene-level tables. It is the boundary between expensive model inference and downstream statistical analysis in the companion repository.

![Downstream aggregation workflow](<../docs/figs/Downstream aggregation workflow: Integrating variant predictions with external biological priors.png>)

## Role in the Thesis Workflow

The thesis computes predicted cis-regulatory genetic variance:

```text
V_G = sum_i 2 p_i (1 - p_i) beta_i^2
```

where `p_i` is allele frequency and `beta_i` is the AlphaGenome-predicted expression effect size. This stage normalizes scored variant tables, preserves RNA track identity, attaches or permutes allele frequencies as required, and aggregates variant-level effects to gene-level `V_G`.

## Main Entrypoint

Use `pipeline_runner.py` for both raw chunk processing and post-processing of existing variant parquet files.

Common environment variables:

- `DATASET_ID`: logical dataset grouping.
- `SAMPLE_ID`: run-specific sample name.
- `ROOT_DIR`: optional experiment root; falls back to scratch or local defaults where configured.

## Core Modules

- `stitcher.py`: reads raw AlphaGenome chunk files.
- `normalizer.py` and `normalisation_helper.py`: standardize variant identifiers and required columns.
- `annotator.py`: attaches allele-frequency information.
- `variant_deduplicator.py`: deduplicates within `gene_id`, `variant_id`, and exact RNA `track_key`.
- `synthetic_variant_downsampler.py`: downsample synthetic/null variants to match observed per-gene counts.
- `permuted_af.py`: assigns observed per-gene allele-frequency distributions to synthetic variants.
- `aggregator.py`: computes gene-level totals, allele-frequency partitions, regional windows, and optional uncertainty summaries.

## Protocol-Aware RNA Identity

AlphaGenome returns multiple RNA tracks. This stage keeps two identifiers explicit:

- `track_key`: stable exact RNA track identity derived from available AlphaGenome metadata.
- `protocol_group`: grouped RNA protocol label, such as `polyA_plus_rna_seq` or `total_rna_seq`.

The same `gene_id` and `variant_id` can legitimately occur in multiple rows when distinct RNA tracks are present. Do not collapse polyA+ and total-RNA tracks unless the downstream analysis explicitly requests that grouping.

## Observed Versus Null Column Contract

Use this contract for thesis-style depletion analyses:

| Dataset type | Total `V_G` | Common/rare partitions | Regional windows |
| --- | --- | --- | --- |
| Observed gnomAD datasets | `vg_predicted` | `vg_common`, `vg_rare` | `vg_promoter_core`, `vg_distal_upstream`, etc. |
| Synthetic/null datasets | `vg_predicted_perm` | `vg_common_perm`, `vg_rare_perm` | `vg_promoter_core_perm`, `vg_distal_upstream_perm`, etc. |

Columns generated only for permutation sanity checks in real observed files are not the biological null model. The biological null comparison uses the separately generated synthetic/null datasets.

## Typical Commands

Observed gnomAD dataset:

```bash
python pipeline_runner.py \
  --variants-parquet /path/to/observed_variant_scores.parquet \
  --variant-out /path/to/observed_variants_dedup.parquet \
  --gene-out /path/to/observed_genes.parquet \
  --deduplicate \
  --permute-af \
  --calc-ci
```

Synthetic/null dataset:

```bash
python pipeline_runner.py \
  --variants-parquet /path/to/synthetic_variant_scores.parquet \
  --real-reference /path/to/observed_variants_dedup.parquet \
  --variant-out /path/to/synthetic_variants_downsampled_perm.parquet \
  --gene-out /path/to/synthetic_genes.parquet \
  --deduplicate \
  --permute-af \
  --calc-ci
```

`--real-reference` is required when a synthetic/null dataset must be downsampled or assigned an observed allele-frequency distribution.

## Outputs

Variant-level outputs contain normalized variant identifiers, allele frequencies, AlphaGenome effect sizes, RNA protocol identity, and optional synthetic/null annotations.

Gene-level outputs contain:

- `protocol_group`
- `n_track_keys`
- `vg_predicted`
- `vg_predicted_perm`
- `vg_common`, `vg_rare`
- regional `V_G` columns relative to the TSS
- optional confidence-interval summaries when `--calc-ci` is enabled
- variant counts and descriptive effect-size summaries

## Interpretation Boundaries

This stage prepares computational predictions for research analysis. It does not establish per-variant pathogenicity. Gene-level `V_G` should be interpreted together with the thesis limitations: modest empirical concordance, unsupported per-variant eQTL concordance in the matched set, and residual sequence-composition confounding in the synthetic null.
