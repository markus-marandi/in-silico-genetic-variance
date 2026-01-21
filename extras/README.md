# extras

Extra scripts for variant generation and analysis.

## generate_synthetic_variants.py

Generate mutation-rate matched synthetic variants.

### workflow

1. **auto-keying**: checks if keyed hail tables exist; if not, automatically runs `load_and_key` on raw tsv files
2. **variant generation**: creates synthetic variants from reference genome within target intervals
3. **mutation rate matching**: matches the 6-class mutation spectrum and methylation distribution of observed variants
4. **cross-dataset disjointness**: ensures generated variants don't overlap between datasets
5. **output generation**: produces hail tables, tsvs, and distribution plots

### usage

```bash
# run with default configuration
python generate_synthetic_variants.py

# or use the example script
./run_synthetic_variants_example.sh

# or configure via environment variables
export INTERMEDIATE_DIR="/path/to/intermediate"
export OUTPUT_DIR="/path/to/output"
python generate_synthetic_variants.py
```

### expected directory structure

under `INTERMEDIATE_DIR` (default: `/mnt/sdb/markus_files/gene_exp/intermediate`):

```
ClinGen_variants.tsv                          # raw variant file
ClinGen_variants.keyed.ht/                    # keyed hail table (auto-created)
ClinGen_gene_set±10kb.tab.bed                 # target intervals

Background_variants.tsv
Background_variants.keyed.ht/
Background_gene_set±10kb.tab.bed
```

### outputs

under `OUTPUT_DIR` (default: `/mnt/sdb/markus_files/gene_exp/51_mutagenesis/ISM_variants`):

```
{dataset}.matched_synthetic.candidates.ht/    # all candidate variants
{dataset}.matched_synthetic.sampled.ht/       # final sampled variants (hail table)
{dataset}.matched_synthetic.sampled.tsv.bgz   # final sampled variants (tsv)
{dataset}.observed_dist.png                   # observed mutation spectrum
{dataset}.sampled_dist.png                    # sampled mutation spectrum
```

### configuration

environment variables with defaults:

- `INTERMEDIATE_DIR`: base directory for inputs (`/mnt/sdb/markus_files/gene_exp/intermediate`)
- `OUTPUT_DIR`: output directory (`/mnt/sdb/markus_files/gene_exp/51_mutagenesis/ISM_variants`)
- `SPARK_CONF_JSON`: spark config json (`/mnt/sdb/markus_files/gene_exp/11_prep_variants/spark_conf.json`)
- `HAIL_HOME`: hail installation (`/mnt/sdb/markus_files/hail`)
- `HAIL_TMP`: temp directory (`/mnt/sdb/markus_files/hail_tmp`)
- `HAIL_LOG`: log file (`/mnt/sdb/markus_files/hail_logs/synthetic_variants.log`)

### integration with path_manager

uses `SyntheticVariantPaths` dataclass from `helpers/path_manager.py` for consistent path management across the project.
