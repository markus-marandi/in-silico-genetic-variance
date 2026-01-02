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


# Background

`python pipeline_runner.py --chunks-dir /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/02_chunks --variants-af /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_variants.tsv --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_gene_set_380.tsv --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/background_genes_20260102.parquet`

`python pipeline_runner.py --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/background_variants_20260102.parquet --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_gene_set_380.tsv --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/background_genes_20260102.parquet`

# Clingen

`python pipeline_runner.py --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/clingen_alphagenome_scores_all_aggs_variantids_long.backfilled.parquet --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/01_inputs/ClinGen_gene_curation_list_GRCh38.ensg.txt --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/clingen_genes_20260102.parquet`

# Background NULL 

`python pipeline_runner.py --chunks-dir /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/02_chunks/chunks --variants-af /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_variants.tsv --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_gene_set_380.tsv `

# Clingen NULL

`python pipeline_runner.py --chunks-dir/cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/02_chunks/chunks --variants-af /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_variants.tsv --gene-list cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/01_inputs/ClinGen_gene_curation_list_GRCh38.ensg.txt`