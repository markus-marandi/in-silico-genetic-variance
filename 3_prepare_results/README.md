single entrypoint
-----------------
- export DATASET_ID, SAMPLE_ID (optionally ROOT_DIR/PDC_TMP).
- run `python pipeline_runner.py --chunks-dir /cfs/.../02_chunks --variants-af /cfs/.../01_inputs/variants.tsv.gz --gnomad-af /cfs/.../gnomad.parquet`.
- optional gene whitelist: `--gene-list /path/to/gene_list.tsv` or place `01_inputs/gene_list.tsv`.

what it does (modules/)
- stitcher.py: lazy scan chunk_*.tsv.gz, normalize gene ids, optional whitelist.
- normalizer.py: ensure CHROM/POS/REF/ALT, gene_tag, method_friendly, variant_id_canonical.
- annotator.py: join AF from initial variants file; optional gnomAD join.
- aggregator.py: gene counts; VG = variance(raw_score) when AF present and non-ISM; promoter/upstream/downstream means; enrich with MANE/GTF/TPM/VGH.
- helpers/external_data_loader.py: loads MANE, GTF, TPM, VGH from $PDC_TMP/initial_data_external.

outputs
- variants parquet: `{sample_id}_variants_{YYYYMMDD}.parquet`
- genes parquet: `{sample_id}_genes_{YYYYMMDD}.parquet`

## To run the pipeline

```bash
parser = argparse.ArgumentParser(description="Stitch, annotate, and summarize chunked scores.")
parser.add_argument("--chunks-dir", type=Path, help="directory containing chunk_*.tsv.gz")
parser.add_argument("--dataset-id", type=str, help="dataset id (preferred)")
parser.add_argument("--sample-id", type=str, help="sample id (preferred)")
parser.add_argument("--root-dir", type=Path, help="override root dir (defaults to $PDC_TMP)")
parser.add_argument("--variant-out", type=Path, help="override variant parquet output path")
parser.add_argument("--gene-out", type=Path, help="override gene parquet output path")
parser.add_argument("--gene-list", type=Path, help="optional gene whitelist (one id per line)")
parser.add_argument("--gnomad-af", type=Path, help="optional gnomAD parquet with CHROM,POS,REF,ALT,AF")
parser.add_argument("--variants-af", type=Path, help="initial variants tsv/tsv.gz with AF (defaults to 01_inputs/variants.tsv.gz)")
parser.add_argument("--variants-parquet", type=Path, help="existing variants parquet to aggregate directly")
```

# Background

`python pipeline_runner.py --chunks-dir /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/02_chunks --variants-af /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_variants.tsv --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_gene_set_380.tsv --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/background_genes_20260102.parquet`

`python pipeline_runner.py --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/background_variants_20260102.parquet --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_gene_set_380.tsv --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/background_genes_20260102.parquet`

# Clingen

`python pipeline_runner.py --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/clingen_alphagenome_scores_all_aggs_variantids_long.backfilled.parquet --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/01_inputs/ClinGen_gene_curation_list_GRCh38.ensg.txt --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/clingen_genes_20260102.parquet`

# Background NULL 

`python pipeline_runner.py --chunks-dir /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/02_chunks/chunks --variants-af /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_variants.tsv --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_gene_set_380.tsv `

`python pipeline_runner.py --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/03_results/dataset5_Background_NULL_variant_level_summary.parquet --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_gene_set_380.tsv --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/03_results/dataset5_Background_NULL_gene_level_summary.parquet`


# Clingen NULL

`python pipeline_runner.py --chunks-dir/cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/02_chunks/chunks --variants-af /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_variants.tsv --gene-list cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/01_inputs/ClinGen_gene_curation_list_GRCh38.ensg.txt`

`python pipeline_runner.py --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/03_results/dataset5_ClinGen_NULL_variant_level_summary.parquet --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/01_inputs/ClinGen_gene_curation_list_GRCh38.ensg.txt --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/03_results/dataset5_ClinGen_NULL_gene_level_summary.parquet`