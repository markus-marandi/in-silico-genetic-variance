#!/bin/bash -l
#SBATCH -A naiss2025-5-479
#SBATCH -J alphagenome_background_NULL_results
#SBATCH -t 01:30:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=96000M
#SBATCH --output=alphagenome_clingen_NULL_results_%j.out
#SBATCH --error=alphagenome_clingen_NULL_results_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=markusmarandi@gmail.com

set -euo pipefail

echo "Job started on $(hostname) at $(date)"

# restore Lmod system defaults
source /opt/cray/pe/cpe/24.11/restore_lmod_system_defaults.sh

ml PDC/24.11
ml miniconda3/25.3.1-1-cpeGNU-24.11
cd /cfs/klemming/home/m/mmarandi/lab/users/mmarandi/alphagenome/in-silico-genetic-variance/3_prepare_results

# use env python directly, avoid activation
ENV_DIR=/cfs/klemming/projects/snic/lappalainen_lab1/users/mmarandi/conda-envs/alphagenome
PYTHON="${ENV_DIR}/bin/python"
if [ ! -x "${PYTHON}" ]; then
  echo "env python not found at ${PYTHON}, using system python" >&2
  PYTHON=python
fi
"${PYTHON}" -V

# ---- install dependency (only if missing) ----
"${PYTHON}" -m pip install --user --upgrade polars
# ---------------------------------------------

# Background NULL
"${PYTHON}" pipeline_runner.py \
  --chunks-dir /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/02_chunks/chunks \
  --gene-list /cfs/klemming/scratch/m/mmarandi/experiments/dataset5/clingen_NULL/01_inputs/ClinGen_gene_curation_list_GRCh38.ensg.txt

echo "Job finished at $(date)"