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
# Unified BAM filter (platform-agnostic). EXPLICIT-ONLY filtering.
# - No default contig whitelist.
# - No default flag filtering.
# - If no filters are provided, it will just sort & index to <in>.filtered.bam.
#
# Usage:
#   bash 03_filter_bam.sh [options] <bam1> [bam2 ...]
#   bash 03_filter_bam.sh [options] -f bam_list.txt
#
# Options:
#   -f FILE            text file with one BAM path per line
#   --contigs FILE     keep ONLY contigs listed (one per line). If omitted, keep all contigs.
#   --exclude-flags N  drop reads with SAM flags N (e.g., 3844). If omitted, keep all flags.
#   --include-flags N  require SAM flags N (e.g., 2 for proper pairs). Optional.
#   --threads N        threads for samtools (default: SLURM_CPUS_PER_TASK | nproc | 8)
#   --suffix STR       output suffix (default: "filtered")
#   --keep-names       preserve original read names (no effect here; placeholder for symmetry)
#   -h | --help        show help
#
# Notes:
# - Output written next to each input as <in>.<suffix>.bam (+ .bai).
# - This script is platform-agnostic; use it for Illumina or HiFi identically.
# - If neither --contigs nor any flag options are provided, the result is a
#   sorted+indexed copy of the input (coordinate order).

THREADS="${SLURM_CPUS_PER_TASK:-$(( $(command -v nproc >/dev/null 2>&1 && nproc) || echo 8 ))}"
LIST=""
CONTIGS=""
EXCL_FLAGS=""
INCL_FLAGS=""
SUFFIX="filtered"

print_help() {
  sed -n '1,200p' "$0" | sed -n '1,80p'
  exit 0
}

ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -f) LIST="${2:?Provide list file}"; shift 2 ;;
    --contigs) CONTIGS="${2:?Provide contigs file}"; shift 2 ;;
    --exclude-flags) EXCL_FLAGS="${2:?Provide integer flag mask}"; shift 2 ;;
    --include-flags) INCL_FLAGS="${2:?Provide integer flag mask}"; shift 2 ;;
    --threads) THREADS="${2:?Provide thread count}"; shift 2 ;;
    --suffix)  SUFFIX="${2:?Provide suffix}"; shift 2 ;;
    --keep-names) shift 1 ;;  # placeholder; no-op
    -h|--help) print_help ;;
    *) ARGS+=("$1"); shift ;;
  esac
done

# Build input list
INPUTS=()
if [[ -n "$LIST" ]]; then
  while IFS= read -r line; do
    [[ -z "$line" ]] && continue
    INPUTS+=("$line")
  done < "$LIST"
fi
INPUTS+=("${ARGS[@]}")
[[ ${#INPUTS[@]} -gt 0 ]] || { echo "[ERR] No BAMs provided."; exit 2; }

# Validate contigs file (optional)
if [[ -n "$CONTIGS" ]]; then
  [[ -s "$CONTIGS" ]] || { echo "[ERR] --contigs file missing or empty: $CONTIGS"; exit 3; }
fi

for INBAM in "${INPUTS[@]}"; do
  [[ -s "$INBAM" ]] || { echo "[ERR] Missing BAM: $INBAM"; exit 4; }
  OUTBAM="${INBAM%.bam}.${SUFFIX}.bam"

  echo "[INFO] Input : $INBAM"
  echo "[INFO] Output: $OUTBAM"
  echo "[INFO] Threads: $THREADS"

  # Ensure input index for sanity checks (idxstats, header)
  if [[ ! -s "${INBAM}.bai" ]]; then
    echo "[INFO] Indexing input BAM ..."
    samtools index -@ "$THREADS" "$INBAM"
  fi

  # Show header contigs & validate requested contigs (if any)
  if [[ -n "$CONTIGS" ]]; then
    echo "[INFO] Validating requested contigs exist in BAM header ..."
    HDR_CONTIGS="$(samtools view -H "$INBAM" | awk -F'\t' '$1=="@SQ"{for(i=1;i<=NF;i++) if($i~"^SN:"){print substr($i,4)}}')"
    MISSING=0
    while read -r C; do
      [[ -z "$C" ]] && continue
      if ! grep -qx "$C" <<< "$HDR_CONTIGS"; then
        echo "[WARN] Contig not present in BAM header: $C"
        MISSING=1
      fi
    done < "$CONTIGS"
    [[ "$MISSING" -eq 1 ]] && echo "[WARN] Some requested contigs are absent; they will be skipped."
  fi

  # Build samtools view command
  # Start with a pass-through; append filters only if explicitly provided.
  VIEW=(samtools view -@ "$THREADS" -b)

  if [[ -n "$EXCL_FLAGS" ]]; then
    VIEW+=(-F "$EXCL_FLAGS")
  fi
  if [[ -n "$INCL_FLAGS" ]]; then
    VIEW+=(-f "$INCL_FLAGS")
  fi

  # Add input BAM
  VIEW+=("$INBAM")

  # Append contig whitelist at the end (if provided)
  if [[ -n "$CONTIGS" ]]; then
    MAPFILE=()
    while IFS= read -r C; do
      [[ -z "$C" ]] && continue
      MAPFILE+=("$C")
    done < "$CONTIGS"
    if [[ ${#MAPFILE[@]} -gt 0 ]]; then
      VIEW+=("${MAPFILE[@]}")
    fi
  fi

  echo "[RUN] ${VIEW[*]} | samtools sort -@ ${THREADS} -o ${OUTBAM} -"
  "${VIEW[@]}" | samtools sort -@ "$THREADS" -o "$OUTBAM" -
  samtools index -@ "$THREADS" "$OUTBAM"

  echo "[STAT] Input idxstats (head):"
  samtools idxstats "$INBAM" | head -n 10
  echo "[STAT] Output idxstats (head):"
  samtools idxstats "$OUTBAM" | head -n 10

  echo "[OK] Wrote: $OUTBAM and ${OUTBAM}.bai"
done

echo "[DONE]"
