# TSS±pad gnomAD v4.1 pipeline

Nextflow pipeline to extract ENSG IDs, build TSS±pad BED from MANE GRCh38, export per-gene gnomAD v4.1 variants with Hail, and merge results.

## Layout
- [`main.nf`](main.nf): pipeline logic and processes.
- [`nextflow.config`](nextflow.config): defaults, profiles, reports.
- Helpers: [`bin/extract_ensg.py`](bin/extract_ensg.py), [`bin/build_tss_bed.py`](bin/build_tss_bed.py), [`bin/hail_export_per_gene.py`](bin/hail_export_per_gene.py), [`bin/merge_variants.py`](bin/merge_variants.py).
- Spark templates: [`conf/spark_local.json`](conf/spark_local.json), [`conf/spark_hpc.json`](conf/spark_hpc.json), [`conf/spark_gcp.json`](conf/spark_gcp.json).
- Runtime example (single config): [`conf/runtime_example.yml`](conf/runtime_example.yml).
- Conda env: [`conf/conda/hail.yml`](conf/conda/hail.yml).

## Inputs
- `--genelist` (required): file containing ENSG IDs (header/no header; any delimiter). ENSG regex: `ENSG\\d+(\\.\\d+)?`; versions are stripped.
- `--mane` (required): MANE GRCh38 summary TSV.
- Optional:
  - `--pad_tss` (default 10000).
  - `--prefix` (default: basename of genelist without .txt/.tsv/.csv/.gz).
  - `--outdir` (default derived; see below).
  - `--gnomad_ht` (default gnomAD v4.1 genomes HT).
  - `--spark_conf`, `--gcs_connector_jar`, `--hail_home`, `--tmp_dir` (default `./tmp`), `--runtime_conf` (YAML blob with the same keys).

## Outdir / prefix derivation
- If `--outdir` is provided, it is used everywhere.
- Else if `DATASET_ID` and `SAMPLE_ID` are set (and optionally `ROOT_DIR`/`PDC_TMP`), output defaults to `$ROOT_DIR/experiments/$DATASET_ID/$SAMPLE_ID/01_inputs/`.
- Else if `--genelist` path contains `/data/initial/`, then `/data/initial/` → `/data/intermediate/` and append `/prefix/`.
- Else default: `./results/<prefix>/`.

## Outputs
- BED: `${outdir}/${prefix}_tss±<pad>.bed`.
- QC: `${outdir}/${prefix}_mapping_qc.tsv`.
- Missing ENSG list: `${outdir}/${prefix}_missing_ensg.txt`.
- Per-gene exports: `${outdir}/per_gene/<gene>.gnomad_v4.1.tsv.bgz`.
- Merge: `${outdir}/${prefix}_variants.tsv.gz`.
- Hail log: `${outdir}/per_gene/hail_export.log`.
- Reports: `${projectDir}/reports/{timeline.html,report.html,trace.txt}`.

## Runtime config (option 2)
Pass `--runtime_conf conf/runtime_example.yml` to set spark_conf, tmp_dir, gnomad_ht, connector jar, hail_home in one file. CLI flags override runtime_conf values.

## Profiles
- `-profile conda`: uses `conf/conda/hail.yml` (includes openjdk17, hail, pyspark, pandas).
- `-profile docker` / `-profile singularity`: use container `hailgenetics/hail:0.2.133`.
- `-profile local` / `hpc` / `gcp`: currently all run locally; pick matching spark_* template and pass via `--spark_conf` or runtime_conf.

## Running examples
- Conda:
  ```
  nextflow run main.nf -profile conda \
    --genelist /path/to/genelist.tsv \
    --mane /path/to/MANE.GRCh38.v1.4.summary.txt \
    --runtime_conf conf/runtime_example.yml
  ```
- Docker:
  ```
  nextflow run main.nf -profile docker \
    --genelist /path/to/genelist.tsv \
    --mane /path/to/MANE.GRCh38.v1.4.summary.txt \
    --spark_conf conf/spark_local.json \
    --gcs_connector_jar /opt/gcs/gcs-connector-hadoop3-latest.jar
  ```

## Notes and validation
- ENSG list must not be empty; pipeline fails fast otherwise.
- MANE symbol is optional; mapping uses longest transcript per ENSG and canonical chromosomes.
- Strand-aware tss: + strand uses transcript start (min coord); - strand uses transcript end (max coord); the bed window centers on that strand-specific tss with ±pad.
- Hail requires a working JVM and gcs connector if accessing `gs://`; if `spark_conf` contains `{GCS_CONNECTOR_JAR}`, you must provide a value or use a container that bundles it.
- Frequency fields are guarded when exporting AF/AC/AN.
