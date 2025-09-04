#!/usr/bin/env bash
set -euo pipefail

# =========================
# EDIT THIS PATH PER SYSTEM
# =========================
# Hard-code BASEDIR here.
BASEDIR="/home/psmc"

# Usage:
#   bash 02_hifi_clean.sh ERR12370312
# Input:
#   BASEDIR/reads/ERR12370312.fastq(.gz)
# Output:
#   BASEDIR/clean_reads/ERR12370312.min1000.q25.p90.hifi.fastq.gz

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <ERR/SRR/DRR>" >&2
    exit 1
fi

RUN="$1"
IN_DIR="${BASEDIR}/reads"
OUT_DIR="${BASEDIR}/clean_reads"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "$OUT_DIR"

# Resolve input file
if [[ -f "${IN_DIR}/${RUN}.fastq.gz" ]]; then
  IN_FASTQ="${IN_DIR}/${RUN}.fastq.gz"
elif [[ -f "${IN_DIR}/${RUN}.fastq" ]]; then
  IN_FASTQ="${IN_DIR}/${RUN}.fastq"
else
  echo "[ERROR] Cannot find ${IN_DIR}/${RUN}.fastq[.gz]" >&2
  exit 2
fi

OUT_FASTQ="${OUT_DIR}/${RUN}.min1000.q25.p90.hifi.fastq.gz"

echo "[INFO] Cleaning HiFi reads"
echo "[INFO] RUN=${RUN}"
echo "[INFO] IN = ${IN_FASTQ}"
echo "[INFO] OUT= ${OUT_FASTQ}"

# Run filtlong with fixed parameters: min length 1000, min Q 25, keep top 90% reads
filtlong --min_length 1000 --min_mean_q 25 --keep_percent 90 "${IN_FASTQ}" \
  | pigz -p "${THREADS}" > "${OUT_FASTQ}"

echo "[OK] Wrote ${OUT_FASTQ}"

