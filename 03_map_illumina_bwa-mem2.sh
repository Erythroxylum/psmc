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
# Usage: bash 03_map_illumina_bwa-mem2.sh <SAMPLE_TAG> <REF_FASTA(.fa/.fna/.fasta[.gz])> <R1.clean.fastq.gz> <R2.clean.fastq.gz> <PLATFORM>
# PLATFORM should be "illumina"
#
# Tips:
# - Submit with plenty of RAM (you said 400G). This script will convert RAM -> speed via large sort buffers.
# - You can override per-thread sort memory with env SORT_MEM_G (default 6G). Example: SORT_MEM_G=8 bash ...
# - Optional downsampling: export TARGET_COV=XX (uses rasusa). Example sbatch: --export=ALL,TARGET_COV=35
# - Uses node-local scratch: $TMPDIR (or /tmp) is cleaned on exit.

# --- args ---
SAMPLE="${1:?SAMPLE_TAG}"
REF="${2:?REF_FASTA}"
R1="${3:?R1.clean.fastq.gz}"
R2="${4:?R2.clean.fastq.gz}"
PLATFORM="${5:?illumina}"

# --- env & paths ---
THREADS="${SLURM_CPUS_PER_TASK:-16}"   # total threads available
MAP_DIR="${BASEDIR}/map"
IDX_DIR="${MAP_DIR}/index"
STATS_FILE="${MAP_DIR}/map_stats.txt"
mkdir -p "${MAP_DIR}" "${IDX_DIR}"

# scratch (node-local) for heavy temps
SCRATCH_ROOT="${TMPDIR:-/tmp}"
TMPD="${SCRATCH_ROOT}/${SAMPLE}.$$"
mkdir -p "${TMPD}"
trap 'rm -rf "$TMPD"' EXIT

# tuning knobs (override via env)
SORT_MEM_G="${SORT_MEM_G:-6}"      # per-thread memory for samtools sort, e.g. 6G (good with 400G RAM)
BWA_BATCH="${BWA_BATCH:-200000000}" # -K batch size for bwa-mem2 (200M ~ good throughput)
ACHIEVED_COV="NA"
CLEAN_DS=0

# tools
for tool in bwa-mem2 samtools; do
  command -v "$tool" >/dev/null 2>&1 || { echo "[ERROR] $tool not found in PATH"; exit 10; }
done

# input checks
[[ -r "${R1}" ]] || { echo "[ERROR] Missing R1: ${R1}" >&2; exit 11; }
[[ -r "${R2}" ]] || { echo "[ERROR] Missing R2: ${R2}" >&2; exit 12; }
[[ -r "${REF}" ]] || { echo "[ERROR] Missing REF: ${REF}" >&2; exit 13; }

# stats header
if [[ ! -f "${STATS_FILE}" ]]; then
  mkdir -p "$(dirname "${STATS_FILE}")"
  echo -e "sample\tplatform\tbam\treads_total\treads_mapped\tpct_mapped\tmean_coverage\ttarget_cov\tachieved_cov" > "${STATS_FILE}"
fi

echo "[INFO] Host: $(hostname)"
echo "[INFO] Threads: ${THREADS}"
echo "[INFO] Scratch: ${TMPD}"
echo "[INFO] SORT_MEM_G: ${SORT_MEM_G}G  (per sort thread)"
echo "[INFO] BWA_BATCH: ${BWA_BATCH} (-K)"

# --- resolve reference into index dir ---
# If REF is gzipped, make a decompressed copy once in IDX_DIR; else symlink.
if [[ "${REF}" =~ \.gz$ ]]; then
  REF_LINK="${IDX_DIR}/${SAMPLE}.ref.fa"
  if [[ ! -f "${REF_LINK}" ]]; then
    echo "[REF] Decompressing ${REF} -> ${REF_LINK}"
    if command -v pigz >/dev/null 2>&1; then pigz -dc "${REF}" > "${REF_LINK}"; else gzip -cd "${REF}" > "${REF_LINK}"; fi
  fi
else
  REF_LINK="${IDX_DIR}/${SAMPLE}.fa"
  [[ -e "${REF_LINK}" ]] || ln -sf "${REF}" "${REF_LINK}"
fi
[[ -e "${REF_LINK}.fai" ]] || samtools faidx "${REF_LINK}"

# bwa-mem2 index (one-time, re-used)
# bwa-mem2 builds multiple *.0123.* files; check one representative.
if ! ls "${REF_LINK}".0123.* >/dev/null 2>&1; then
  echo "[IDX] bwa-mem2 index ${REF_LINK}"
  bwa-mem2 index "${REF_LINK}"
fi

# --- optional downsample with rasusa to TARGET_COV ---
if [[ -n "${TARGET_COV:-}" ]]; then
  command -v rasusa >/dev/null 2>&1 || { echo "[ERROR] rasusa not found but TARGET_COV set"; exit 14; }
  echo "[DS] TARGET_COV=${TARGET_COV}x → computing genome size from ${REF_LINK}.fai"
  GENOME_SIZE=$(awk '{s+=$2} END{print s}' "${REF_LINK}.fai")
  DS_DIR="${MAP_DIR}/downsampled/${SAMPLE}"
  mkdir -p "${DS_DIR}"
  DS_PREFIX="${DS_DIR}/${SAMPLE}_cov${TARGET_COV}"

  echo "[DS] rasusa → ${DS_PREFIX}_{1,2}.fastq.gz"
  rasusa reads "${R1}" "${R2}" -g "${GENOME_SIZE}" -c "${TARGET_COV}" \
    -o "${DS_PREFIX}_1.fastq.gz" -o "${DS_PREFIX}_2.fastq.gz" -s 42

  R1="${DS_PREFIX}_1.fastq.gz"
  R2="${DS_PREFIX}_2.fastq.gz"
  CLEAN_DS=1

  bases=$(zcat "${DS_PREFIX}_1.fastq.gz" "${DS_PREFIX}_2.fastq.gz" | awk 'NR%4==2 {sum+=length($0)} END{print sum+0}')
  if [[ -n "${GENOME_SIZE}" && "${GENOME_SIZE}" -gt 0 ]]; then
    ACHIEVED_COV=$(awk -v b="$bases" -v g="$GENOME_SIZE" 'BEGIN{printf "%.2f", b/g}')
    echo "[DS] Achieved coverage: target=${TARGET_COV}x, actual=${ACHIEVED_COV}x"
  fi
fi

# --- thread split ---
# Give most threads to bwa-mem2; keep enough for sort/markdup parallelism.
BWA_T=$(( THREADS >= 24 ? 16 : (THREADS*2/3) ))
SRT_T=$(( THREADS - BWA_T ))
if (( SRT_T < 6 )); then SRT_T=6; BWA_T=$(( THREADS - SRT_T )); fi
echo "[TUNE] bwa-mem2 threads=${BWA_T}, sort/markdup threads=${SRT_T}"

# sort memory per thread (e.g., 6G) → total sort buffers ≈ SRT_T * SORT_MEM_G
SORT_MEM="${SORT_MEM_G}G"
echo "[TUNE] samtools sort per-thread memory=${SORT_MEM} (total ~ $((SRT_T * SORT_MEM_G))G)"

# --- outputs ---
POS_BAM="${TMPD}/${SAMPLE}.ill.pos.bam"                      # intermediate on scratch (uncompressed for speed)
MARKDUP_BAM="${TMPD}/${SAMPLE}.ill.sorted.markdup.bam"       # final on scratch, then copied out
FINAL_OUT="${MAP_DIR}/${SAMPLE}.ill.sorted.markdup.bam"

# --- mapping pipeline ---
echo "[MAP] bwa-mem2 mem -K ${BWA_BATCH} -t ${BWA_T}  | fixmate | sort (t=${SRT_T}, -m ${SORT_MEM}, -T ${TMPD})"
bwa-mem2 mem -K "${BWA_BATCH}" -t "${BWA_T}" \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
  "${REF_LINK}" "${R1}" "${R2}" \
| samtools fixmate -@ "${SRT_T}" -m - - \
| samtools sort   -@ "${SRT_T}" -m "${SORT_MEM}" -l 0 -T "${TMPD}/sort" -o "${POS_BAM}" -

echo "[DUP] markdup (t=${SRT_T}) -> ${MARKDUP_BAM}"
samtools markdup -@ "${SRT_T}" -T "${TMPD}/markdup" -r -s "${POS_BAM}" "${MARKDUP_BAM}"
rm -f "${POS_BAM}"

samtools index -@ "${SRT_T}" "${MARKDUP_BAM}"

# copy final BAM + index back to MAP_DIR
cp -f "${MARKDUP_BAM}" "${FINAL_OUT}"
cp -f "${MARKDUP_BAM}.bai" "${FINAL_OUT}.bai"
echo "[OK] ${FINAL_OUT}"

# --- QC stats (robust parsing) ---
flagstat=$(samtools flagstat -@ "${SRT_T}" "${FINAL_OUT}")
reads_total=$(echo "$flagstat"   | awk '/in total/ {print $1; exit}')
reads_mapped=$(echo "$flagstat"  | awk '/ mapped \(/ {print $1; exit}')
mapped_pct=$(echo "$flagstat"    | awk -F'[()%]' '/ mapped \(/ {gsub(/^ +| +$/,"",$2); print $2; exit}')
mean_cov=$(samtools depth -@ "${SRT_T}" -a "${FINAL_OUT}" | awk '{sum+=$3} END{ if(NR>0) printf "%.2f", sum/NR; else print "0.00"}')

# refresh stats row
tmp="${STATS_FILE}.tmp.$$"
awk -v s="${SAMPLE}" -F'\t' 'NR==1{print;next} $1!=s' "${STATS_FILE}" > "${tmp}" || true
echo -e "${SAMPLE}\t${PLATFORM}\t${FINAL_OUT}\t${reads_total}\t${reads_mapped}\t${mapped_pct}\t${mean_cov}\t${TARGET_COV:-NA}\t${ACHIEVED_COV:-NA}" >> "${tmp}"
mv "${tmp}" "${STATS_FILE}"

# cleanup downsample unless KEEP_DOWNSAMPLED=1
if [[ "${CLEAN_DS:-0}" -eq 1 && "${KEEP_DOWNSAMPLED:-0}" -ne 1 ]]; then
  rm -f "${DS_PREFIX}_1.fastq.gz" "${DS_PREFIX}_2.fastq.gz" 2>/dev/null || true
fi

echo "[DONE] Mapping complete."


