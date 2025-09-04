#!/usr/bin/env bash
set -euo pipefail

# =========================
# EDIT THIS PATH PER SYSTEM
# =========================
# Hard-code BASEDIR here.
BASEDIR="/home/psmc"

# Usage: bash 01_fetch_reads_sra.sh SRR5313982

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <SRA_RUN_ID>" >&2
    exit 1
fi

RUN_ID="$1"
OUTDIR="${BASEDIR}/reads"
SRACACHE="${BASEDIR}/reads/sra_cache"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "$OUTDIR" "$SRACACHE"

echo "[INFO] Fetching $RUN_ID"
echo "[INFO] Output directory: $OUTDIR"
echo "[INFO] Threads: $THREADS"
echo "[INFO] SRA cache: $SRACACHE"

# Point SRA cache to scratch
vdb-config -s "/repository/user/main/public/root=${SRACACHE}" || true

# Put tmp on scratch so no big files land in $PWD/$HOME
export TMPDIR="${SRACACHE}/tmp"
mkdir -p "$TMPDIR"

prefetch "$RUN_ID" -X 36G || echo "[WARN] prefetch failed or already cached, continuing…"

# Use --temp if available, otherwise rely on TMPDIR
if fasterq-dump --help 2>&1 | grep -q -- '--temp'; then
    fasterq-dump "$RUN_ID" --split-files -O "$OUTDIR" -e "$THREADS" --temp "$TMPDIR"
else
    fasterq-dump "$RUN_ID" --split-files -O "$OUTDIR" -e "$THREADS"
fi

echo "[INFO] Download complete:"
ls -lh "$OUTDIR" | grep "$RUN_ID" || echo "[WARN] No files matched ${RUN_ID} in ${OUTDIR}"

