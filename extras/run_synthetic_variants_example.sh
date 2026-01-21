#!/bin/bash
# example script to run unified synthetic variant generation

# configure paths via environment variables
export INTERMEDIATE_DIR="/mnt/sdb/markus_files/gene_exp/intermediate"
export OUTPUT_DIR="/mnt/sdb/markus_files/gene_exp/51_mutagenesis/ISM_variants"
export SPARK_CONF_JSON="/mnt/sdb/markus_files/gene_exp/11_prep_variants/spark_conf.json"
export HAIL_HOME="/mnt/sdb/markus_files/hail"
export HAIL_TMP="/mnt/sdb/markus_files/hail_tmp"
export HAIL_LOG="/mnt/sdb/markus_files/hail_logs/synthetic_variants.log"

# run the unified workflow
# will automatically:
# 1. check if keyed tables exist
# 2. run load_and_key if raw tsvs exist but keyed tables don't
# 3. generate synthetic variants with mutation rate matching
python generate_synthetic_variants.py
