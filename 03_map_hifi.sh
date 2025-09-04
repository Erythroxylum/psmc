#!/usr/bin/env bash
set -euo pipefail
#
# =========================
# EDIT THIS PATH PER SYSTEM
# =========================
# Hard-code BASEDIR here.
BASEDIR="/home/psmc"
#
#
# Usage: bash 03_map_hifi.sh <SAMPLE_TAG> <REF_FASTA(.fa/.fna/.fasta[.gz])> <HiFi.clean.fastq.gz> <PLATFORM>
# PLATFORM should be "hifi"
#
# Optional env:
#   TARGET_COV=20         # downsample to ~20× via rasusa (uses reference size)
#   KEEP_DOWNSAMPLED=1    # keep rasusa outputs; default is to remove after mapping

# --- args ---
SAMPLE="${1:?SAMPLE_TAG}"
REF="${2:?REF_FASTA}"
READS_IN="${3:?HiFi.clean.fastq.gz}"
PLATFORM="${4:?hifi}"

# --- env & paths ---
THREADS="${SLURM_CPUS_PER_TASK:-8}"
MAP_DIR="${BASEDIR}/map"
IDX_DIR="${MAP_DIR}/index"
TMP_DIR="${MAP_DIR}/tmp/${SAMPLE}"
STATS_FILE="${MAP_DIR}/map_stats.txt"
mkdir -p "${MAP_DIR}" "${IDX_DIR}" "${TMP_DIR}"

# Ensure stats file exists with header
if [[ ! -f "${STATS_FILE}" ]]; then
  mkdir -p "$(dirname "${STATS_FILE}")"
  echo -e "sample\tplatform\tbam\treads_total\treads_mapped\tpct_mapped\tmean_coverage\ttarget_cov\tachieved_cov" > "${STATS_FILE}"
fi

# --- tool checks ---
for tool in minimap2 samtools; do
  command -v "$tool" >/dev/null 2>&1 || { echo "[ERROR] $tool not found in PATH"; exit 10; }
done
command -v rasusa >/dev/null 2>&1 || true
command -v pigz >/dev/null 2>&1 || true

# --- input checks ---
[[ -r "${READS_IN}" ]] || { echo "[ERROR] Missing/unreadable reads: ${READS_IN}" >&2; exit 11; }
[[ -r "${REF}"     ]] || { echo "[ERROR] Missing/unreadable REF: ${REF}"     >&2; exit 12; }

# --- helper: fix malformed FASTQ (seq/qual length mismatch) ---
fix_fastq () {
  local in="$1" out="$2" thr="${3:-4}"
  echo "[fix-fastq] Filtering malformed records: ${in} → ${out}"
  local DECOMP COMP
  if command -v pigz >/dev/null 2>&1; then
    DECOMP="pigz -dc"; COMP="pigz -p ${thr}"
  else
    DECOMP="gzip -cd"; COMP="gzip -c"
  fi
  $DECOMP "$in" \
  | awk '
      NR%4==1 {h=$0}
      NR%4==2 {s=$0}
      NR%4==3 {p=$0}
      NR%4==0 {
        q=$0
        if (length(s)==length(q)) {
          print h; print s; print p; print q
        } else {
          bad++
        }
      }
      END {
        if (bad>0) {
          printf("[fix-fastq] Dropped %d bad records\n", bad) > "/dev/stderr"
        }
      }
    ' \
  | $COMP > "$out"
  echo "[fix-fastq] Wrote: $out"
}

# --- resolve reference into index dir ---
# If REF is gzipped, decompress once to IDX_DIR; else symlink.
if [[ "${REF}" =~ \.gz$ ]]; then
  REF_LINK="${IDX_DIR}/${SAMPLE}.ref.fa"
  if [[ ! -f "${REF_LINK}" ]]; then
    echo "[REF] Decompressing ${REF} -> ${REF_LINK}"
    if command -v pigz >/dev/null 2>&1; then
      pigz -dc "${REF}" > "${REF_LINK}"
    else
      gzip -cd "${REF}" > "${REF_LINK}"
    fi
  fi
else
  REF_LINK="${IDX_DIR}/${SAMPLE}.fa"
  [[ -e "${REF_LINK}" ]] || ln -sf "${REF}" "${REF_LINK}"
fi

# --- index reference ---
[[ -e "${REF_LINK}.fai" ]] || samtools faidx "${REF_LINK}"

# --- working read path (may be replaced) ---
READS="${READS_IN}"
CLEAN_DS=0
FIXED_MADE=0
ACHIEVED_COV="NA"

# --- optional downsample with rasusa to TARGET_COV ---
if [[ -n "${TARGET_COV:-}" ]]; then
  command -v rasusa >/dev/null 2>&1 || { echo "[ERROR] rasusa not found but TARGET_COV set"; exit 13; }

  echo "[DS] TARGET_COV=${TARGET_COV}x → computing genome size from ${REF_LINK}"
  [[ -e "${REF_LINK}.fai" ]] || samtools faidx "${REF_LINK}"
  GENOME_SIZE=$(awk '{s+=$2} END{print s}' "${REF_LINK}.fai")
  echo "[DS] genome_size=${GENOME_SIZE}"

  DS_DIR="${MAP_DIR}/downsampled/${SAMPLE}"
  mkdir -p "${DS_DIR}"
  DS_PREFIX="${DS_DIR}/${SAMPLE}_cov${TARGET_COV}"
  DS_OUT="${DS_PREFIX}.fastq.gz"

  run_rasusa () {
    echo "[DS] rasusa → ${DS_OUT}"
    rasusa reads \
      "${READS}" \
      -g "${GENOME_SIZE}" \
      -c "${TARGET_COV}" \
      -o "${DS_OUT}" \
      -s 42
  }

  # Try rasusa; if it fails with parse error, fix and retry once
  if ! run_rasusa 2> "${TMP_DIR}/rasusa.err"; then
    if grep -qiE 'Failed to parse record|unable to gather read lengths|Sequence length.*quality length' "${TMP_DIR}/rasusa.err"; then
      echo "[DS] rasusa failed due to malformed FASTQ. Auto-fixing and retrying once…"
      FIXED_FASTQ="${TMP_DIR}/$(basename "${READS_IN%.gz}").fixed.fastq.gz"
      fix_fastq "${READS_IN}" "${FIXED_FASTQ}" "${THREADS}"
      READS="${FIXED_FASTQ}"
      FIXED_MADE=1
      run_rasusa
    else
      echo "[DS][WARN] rasusa failed for another reason; proceeding WITHOUT downsampling."
      DS_OUT=""
    fi
  fi

  # If downsampled file exists, swap inputs and compute achieved coverage
  if [[ -n "${DS_OUT}" && -s "${DS_OUT}" ]]; then
    READS="${DS_OUT}"
    CLEAN_DS=1
    bases=$(zcat "${READS}" | awk 'NR%4==2 {sum+=length($0)} END{print sum+0}')
    if [[ "${GENOME_SIZE}" -gt 0 ]]; then
      ACHIEVED_COV=$(awk -v b="$bases" -v g="$GENOME_SIZE" 'BEGIN{printf "%.2f", b/g}')
      echo "[DS] Achieved coverage: target=${TARGET_COV}x, actual=${ACHIEVED_COV}x"
    fi
  fi
fi

# --- outputs ---
BAM="${MAP_DIR}/${SAMPLE}.hifi.sorted.bam"

echo "[MAP] minimap2 -ax map-hifi | sort -> ${BAM}"
minimap2 -t "${THREADS}" -ax map-hifi "${REF_LINK}" "${READS}" \
  | samtools sort -@ "${THREADS}" -o "${BAM}" -
samtools index -@ "${THREADS}" "${BAM}"

# --- QC stats (robust parsing) ---
flagstat=$(samtools flagstat -@ "${THREADS}" "${BAM}")
reads_total=$(echo "$flagstat" | awk '/in total/ {print $1; exit}')
reads_mapped=$(echo "$flagstat" | awk '/ mapped \(/ {print $1; exit}')
mapped_pct=$(echo "$flagstat" | awk -F'[()%]' '/ mapped \(/ {gsub(/^ +| +$/,"",$2); print $2; exit}')
mean_cov=$(samtools depth -@ "${THREADS}" -a "${BAM}" | awk '{sum+=$3} END{ if(NR>0) printf "%.2f", sum/NR; else print "0.00"}')

# --- write/refresh stats row ---
mkdir -p "$(dirname "${STATS_FILE}")"
if [[ ! -f "${STATS_FILE}" ]]; then
  echo -e "sample\tplatform\tbam\treads_total\treads_mapped\tpct_mapped\tmean_coverage\ttarget_cov\tachieved_cov" > "${STATS_FILE}"
fi
tmp="${STATS_FILE}.tmp.$$"
awk -v s="${SAMPLE}" -F'\t' 'NR==1{print;next} $1!=s' "${STATS_FILE}" > "${tmp}" || true
echo -e "${SAMPLE}\t${PLATFORM}\t${BAM}\t${reads_total}\t${reads_mapped}\t${mapped_pct}\t${mean_cov}\t${TARGET_COV:-NA}\t${ACHIEVED_COV:-NA}" >> "${tmp}"
mv "${tmp}" "${STATS_FILE}"

echo "[OK] ${BAM}"

# --- cleanup temp (downsampled & fixed) unless KEEP_DOWNSAMPLED=1 ---
if [[ "${KEEP_DOWNSAMPLED:-0}" -ne 1 ]]; then
  [[ -n "${DS_OUT:-}" ]] && rm -f "${DS_OUT}" 2>/dev/null || true
  if [[ "${FIXED_MADE}" -eq 1 ]]; then
    rm -f "${FIXED_FASTQ}" 2>/dev/null || true
  fi
fi

