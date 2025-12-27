# 1. Allocate resources
salloc \
  -A naiss2025-5-479 \
  -J alphagenome_filtering_int \
  -t 1:00:00 \
  -p shared \
  -n 1 \
  -c 2 \
  --mem=64000M

# when the prompt changes to the compute node, start an interactive shell
srun --pty bash -l

# restore Lmod defaults
source /opt/cray/pe/cpe/24.11/restore_lmod_system_defaults.sh
ml PDC/24.11
ml miniconda3/25.3.1-1-cpeGNU-24.11

cd /cfs/klemming/home/m/mmarandi/lab/users/mmarandi/alphagenome

# env vars similar to your sbatch --export
export ENV_FILE=/cfs/klemming/home/m/mmarandi/lab/users/mmarandi/alphagenome/src/.env
export API_KEY_VAR=API_KEY_PERSONAL
export JOB_INDEX=0
export JOB_TOTAL=3

ENV_DIR=/cfs/klemming/projects/snic/lappalainen_lab1/users/mmarandi/conda-envs/alphagenome
PYTHON="${ENV_DIR}/bin/python"

PYTHONUNBUFFERED=1 "$PYTHON" src/31_alphagenome_batch_scoring_batch_window.py