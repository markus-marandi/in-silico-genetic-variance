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

`python pipeline_runner.py --chunks-dir /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/02_chunks --variants-af /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_variants.tsv --gnomad-af /cfs/.../gnomad.parquet`