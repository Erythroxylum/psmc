#!/bin/bash
#SBATCH -p sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o log.%j.out
#SBATCH -e log.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<EMAIL>
#
#
# =========================
# EDIT SCRIPTDIR PATH PER SYSTEM
# =========================
# Hard-code SCRIPTDIR here.
SCRIPTDIR="$HOME/scripts/psmc/03_map"
#
# Usage:
#   sbatch --job-name=map-ILL-<TAG> step03_map_illumina_bwa-mem2.sbatch <TAG> <REF> <R1.clean.gz> <R2.clean.gz> <PLATFORM>
# Example (with env knobs passed through): 
#   sbatch --export=ALL,TARGET_COV=35,SORT_MEM_G=8,BWA_BATCH=300000000 \
#     --job-name=map-ILL-Pfurb step03_map_illumina_bwa-mem2.sbatch Pfurb ref.fa R1.fq.gz R2.fq.gz illumina

set -euo pipefail

[[ $# -eq 5 ]] || { echo "Usage: $0 <TAG> <REF> <R1.clean.gz> <R2.clean.gz> <PLATFORM>"; exit 2; }
TAG="$1"; REF="$2"; R1="$3"; R2="$4"; PLATFORM="$5"

SCRIPT="$SCRIPTDIR/03_map_illumina_bwa-mem2.sh"

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate psmc

echo "[INFO] Host: $(hostname)"
echo "[INFO] TMPDIR: ${TMPDIR:-/tmp}"
echo "[INFO] CPUs: ${SLURM_CPUS_PER_TASK}"
echo "[INFO] MEM: ${SLURM_MEM_PER_NODE:-${SLURM_MEM_PER_CPU:-unknown}}"
echo "[INFO] TARGET_COV=${TARGET_COV:-unset}  SORT_MEM_G=${SORT_MEM_G:-6}  BWA_BATCH=${BWA_BATCH:-200000000}"

bash "$SCRIPT" "$TAG" "$REF" "$R1" "$R2" "$PLATFORM"

