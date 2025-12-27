#!/bin/bash -l
#SBATCH -A naiss2025-5-479
#SBATCH -J alphagenome_window
#SBATCH -t 60:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=14000M
#SBATCH --output=alphagenome_window_%j.out
#SBATCH --error=alphagenome_window_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=markusmarandi@gmail.com

echo "Job started on $(hostname) at $(date)"

# restore Lmod system defaults
source /opt/cray/pe/cpe/24.11/restore_lmod_system_defaults.sh

ml PDC/24.11
ml miniconda3/25.3.1-1-cpeGNU-24.11
cd /cfs/klemming/home/m/mmarandi/lab/users/mmarandi/alphagenome/master-thesis-pipeline/2_score_variants_with_alphagenome

# use env python directly, avoid activation
ENV_DIR=/cfs/klemming/projects/snic/lappalainen_lab1/users/mmarandi/conda-envs/alphagenome
PYTHON="${ENV_DIR}/bin/python"
if [ ! -x "${PYTHON}" ]; then
echo "env python not found at ${PYTHON}, using system python" >&2
PYTHON=python
fi
"${PYTHON}" -V

# load .env and normalize api key var name
ENV_FILE="${ENV_FILE:-.env}"
if [ ! -f "${ENV_FILE}" ] && [ -f src/.env ]; then
ENV_FILE="src/.env"
fi
if [ -f "${ENV_FILE}" ]; then
set -a
. "${ENV_FILE}"
set +a
else
echo "warning: env file ${ENV_FILE} not found" >&2
fi
API_KEY_VAR="${API_KEY_VAR:-API_KEY_PERSONAL}"
if [ -n "${API_KEY_OVERRIDE:-}" ]; then
SELECTED_API_KEY="${API_KEY_OVERRIDE}"
else
SELECTED_API_KEY="${!API_KEY_VAR:-}"
fi
if [ -z "${SELECTED_API_KEY:-}" ]; then
echo "error: api key not found via ${API_KEY_VAR} or API_KEY_OVERRIDE" >&2
exit 1
fi
export API_KEY_PERSONAL="${SELECTED_API_KEY}"

# require dataset/sample to drive dataset-centric layout
if [ -z "${DATASET_ID:-}" ] || [ -z "${SAMPLE_ID:-}" ]; then
  echo "error: DATASET_ID and SAMPLE_ID must be set" >&2
  exit 1
fi
export DATASET_ID
export SAMPLE_ID

# no package upgrades here; assume env is pre-provisioned

# default job fan-out; override via --export
: "${JOB_INDEX:=0}"
: "${JOB_TOTAL:=1}"
export JOB_INDEX
export JOB_TOTAL
echo "using api key var ${API_KEY_VAR}, job $((JOB_INDEX + 1))/${JOB_TOTAL}"

PYTHONUNBUFFERED=1 "${PYTHON}" /cfs/klemming/home/m/mmarandi/lab/users/mmarandi/alphagenome/master-thesis-pipeline/2_score_variants_with_alphagenome/alphagenome_batch_scoring_batch_window.py

status=$?
echo "Job finished with status ${status} at $(date)"
exit ${status}

# TO RUN IT
# ENV_PATH=/cfs/klemming/home/m/mmarandi/lab/users/mmarandi/alphagenome/master-thesis-pipeline/2_score_variants_with_alphagenome/.env
# sbatch --export=ENV_FILE=${ENV_PATH},API_KEY_VAR=API_KEY_PERSONAL,JOB_INDEX=0,JOB_TOTAL=3,DATASET_ID=dataset4,SAMPLE_ID=background  run_alphagenome_batch_window.sh
# sbatch --export=ENV_FILE=${ENV_PATH},API_KEY_VAR=API_KEY_1674_PERSONAL,JOB_INDEX=1,JOB_TOTAL=3,DATASET_ID=dataset4,SAMPLE_ID=background run_alphagenome_batch_window.sh
# sbatch --export=ENV_FILE=${ENV_PATH},API_KEY_VAR=API_KEY_SCILIFELAB,JOB_INDEX=2,JOB_TOTAL=3,DATASET_ID=dataset4,SAMPLE_ID=background run_alphagenome_batch_window.sh