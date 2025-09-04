#!/usr/bin/env bash
set -euo pipefail

# =========================
# EDIT THIS PATH PER SYSTEM
# =========================
# Hard-code BASEDIR here.
BASEDIR="/home/psmc"

# Usage: bash fastp.sh SRR123456

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <SRR_RUN_ID>" >&2
    exit 1
fi

RUN_ID="$1"
READS_DIR="${BASEDIR}/reads"
CLEAN_DIR="${BASEDIR}/clean_reads"

R1="${READS_DIR}/${RUN_ID}_1.fastq"
R2="${READS_DIR}/${RUN_ID}_2.fastq"

# Check input files
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "[ERROR] One or both input files not found:"
    echo "  $R1"
    echo "  $R2"
    exit 2
fi

# Make output dir
mkdir -p "$CLEAN_DIR"

# Outputs
OUT_R1="${CLEAN_DIR}/${RUN_ID}_R1.clean.fastq.gz"
OUT_R2="${CLEAN_DIR}/${RUN_ID}_R2.clean.fastq.gz"
HTML_REPORT="${CLEAN_DIR}/${RUN_ID}.fastp.html"
JSON_REPORT="${CLEAN_DIR}/${RUN_ID}.fastp.json"

echo "[INFO] Running fastp on ${RUN_ID}"
fastp \
    -i "$R1" \
    -I "$R2" \
    -o "$OUT_R1" \
    -O "$OUT_R2" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --trim_front1=5 \
    --trim_front2=5 \
    --cut_mean_quality=20 \
    --cut_front \
    --cut_tail \
    --trim_poly_g \
    --trim_poly_x \
    --length_required=35 \
    --dedup \
    --thread "${SLURM_CPUS_PER_TASK:-1}" \
    --html "$HTML_REPORT" \
    --json "$JSON_REPORT"

echo "[INFO] fastp complete for ${RUN_ID}"
echo "[INFO] Cleaned reads in: $CLEAN_DIR"

