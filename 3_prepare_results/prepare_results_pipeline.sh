#!/bin/bash -l
#SBATCH -A naiss2025-5-479
#SBATCH -J alphagenome_results_pipeline
#SBATCH -t 03:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=128000M
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=markusmarandi@gmail.com

set -euo pipefail

echo "job started on $(hostname) at $(date)"
echo "sbatch job id: $SLURM_JOB_ID"

module reset
ml PDC/24.11
ml miniconda3/25.3.1-1-cpeGNU-24.11

# set paths
WORK_DIR=/cfs/klemming/home/m/mmarandi/lab/users/mmarandi/alphagenome/in-silico-genetic-variance/3_prepare_results
cd "$WORK_DIR"

ENV_DIR=/cfs/klemming/projects/snic/lappalainen_lab1/users/mmarandi/conda-envs/alphagenome
PYTHON="${ENV_DIR}/bin/python"

# create logs dir if needed
mkdir -p logs

# verify python exists
if [ ! -f "$PYTHON" ]; then
  echo "error: python not found at $PYTHON"
  exit 1
fi

echo "starting pipeline runs at $(date)"
echo ""

# dataset 3: clingen (real gnomad variants)
echo "===== dataset3 clingen (real) ====="
echo "start: $(date)"
$PYTHON pipeline_runner.py \
  --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/clingen_alphagenome_scores_all_aggs_variantids_long.backfilled.parquet \
  --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/01_inputs/ClinGen_gene_curation_list_GRCh38.ensg.txt \
  --variant-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/ClinGen_HI_Gnomad_variants_dedup_14022026.parquet \
  --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/ClinGen_HI_Gnomad_genes_14022026.parquet \
  --deduplicate \
  --permute-af \
  --calc-ci

echo "done: $(date)"
echo ""

# dataset 4: background (real gnomad variants)
echo "===== dataset4 background (real) ====="
echo "start: $(date)"
$PYTHON pipeline_runner.py \
  --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/background_variants_20260102.parquet \
  --variant-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/Background_Gnomad_variants_dedup_perm_14022026.parquet \
  --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/Background_Gnomad_genes_14022026.parquet \
  --deduplicate \
  --permute-af \
  --calc-ci

echo "done: $(date)"
echo ""

# dataset 5: background null (synthetic variants)
echo "===== dataset5 background_NULL (synthetic) ====="
echo "start: $(date)"
$PYTHON pipeline_runner.py \
  --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/03_results/dataset5_Background_NULL_variant_level_summary.parquet \
  --real-reference /cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/03_results/Background_Gnomad_variants_dedup_perm_14022026.parquet \
  --variant-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/03_results/Background_Synth_variants_downsampled_perm_14022026.parquet \
  --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/background_NULL/03_results/Background_Synth_genes_14022026.parquet \
  --deduplicate \
  --permute-af \
  --calc-ci

echo "done: $(date)"
echo ""

# dataset 5: clingen null (synthetic variants)
echo "===== dataset5 clingen_NULL (synthetic) ====="
echo "start: $(date)"
$PYTHON pipeline_runner.py \
  --variants-parquet /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/03_results/dataset5_ClinGen_NULL_variant_level_summary.parquet \
  --real-reference /cfs/klemming/scratch/m/mmarandi/experiments/dataset3/clingen/03_results/ClinGen_HI_Gnomad_variants_dedup_14022026.parquet \
  --variant-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/03_results/ClinGen_HI_Synth_variants_downsampled_perm_14022026.parquet \
  --gene-out /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/03_results/ClinGen_HI_Synth_genes_14022026.parquet \
  --deduplicate \
  --permute-af \
  --calc-ci

echo "done: $(date)"
echo ""

echo "===== all pipelines completed successfully ====="
echo "end: $(date)"