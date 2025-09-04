#!/usr/bin/env bash
set -euo pipefail

# 06_mask_coverage_outliers.sh
#
# Generate BED masks of sites with *outlier* coverage for a list of BAM files.
# These masks can be fed to `bcftools consensus -m` to set those regions to N.
#
# USAGE
#   # Minimal (defaults: --min 10 --max 60)
#   bash 06_mask_coverage_outliers.sh -f bam_list.txt
#
#   # Custom coverage and quality filters
#   bash 06_mask_coverage_outliers.sh -f bam_list.txt --min 8 --max 80 -q 1 -Q 13
#
#   # Faster with threads and custom output dir
#   bash 06_mask_coverage_outliers.sh -f bam_list.txt -t 8 -o coverage_masks
#
# INPUT
#   -f, --file LIST      Plain-text file of BAM paths (one per line; absolute or relative).
#
# COVERAGE OUTLIER DEFINITION
#   --min INT            Minimum depth to KEEP (mask depth < MIN).   [default: 10]
#   --max INT            Maximum depth to KEEP (mask depth > MAX).   [default: 60]
#
# DEPTH QUALITY FILTERS (applied before outlier detection)
#   -q INT               samtools depth: minimum mapping quality.    [default: 0]
#   -Q INT               samtools depth: minimum base quality.       [default: 0]
#
# PERFORMANCE / OUTPUT
#   -t, --threads INT    Threads for samtools index/depth.           [default: 1]
#   -o DIR               Output directory for BED masks.             [default: ./coverage_masks]
#
# OUTPUT
#   For each input BAM: <OUTDIR>/<sample>.mask.bed
#   where <sample> is the BAM basename without the .bam extension.
#   The BED reports *regions to MASK* (depth < MIN OR depth > MAX), merged.
#
# NOTES
#   * Requires: samtools, bedtools.
#   * The BAM must be coordinate-sorted; the script will create a .bai index if missing.
#   * We sort sites before `bedtools merge` to satisfy BED sorting expectations.
#   * Positions with zero coverage are included via `samtools depth -a`.
#   * Example downstream use:
#       bcftools consensus -f ref.fa -H I -s SAMPLE -m sample.mask.bed sample.vcf.gz > sample.iupac.fa
#
# HELP
#   -h, --help           Print this usage block.

# -------------------- arg parse --------------------

LIST=""
MIN=10
MAX=60
OUTDIR="./coverage_masks"
MAPQ=0
BASEQ=0
THREADS=1

print_help() { sed -n '4,49p' "$0" | sed 's/^# \{0,1\}//'; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -f|--file)   LIST="${2:?}"; shift 2 ;;
    --min)       MIN="${2:?}";  shift 2 ;;
    --max)       MAX="${2:?}";  shift 2 ;;
    -o)          OUTDIR="${2:?}"; shift 2 ;;
    -q)          MAPQ="${2:?}"; shift 2 ;;
    -Q)          BASEQ="${2:?}"; shift 2 ;;
    -t|--threads) THREADS="${2:?}"; shift 2 ;;
    -h|--help)   print_help; exit 0 ;;
    *) echo "[ERR] Unknown argument: $1" >&2; exit 2 ;;
  esac
done

# -------------------- checks --------------------

need() { command -v "$1" >/dev/null 2>&1 || { echo "[ERR] Missing tool: $1" >&2; exit 3; }; }
need samtools
need bedtools

[[ -n "$LIST" ]] || { echo "[ERR] Provide -f/--file bam_list.txt" >&2; exit 2; }
[[ -s "$LIST" ]] || { echo "[ERR] BAM list file is empty: $LIST" >&2; exit 2; }

# integers
[[ "$MIN" =~ ^[0-9]+$ ]] || { echo "[ERR] --min must be integer" >&2; exit 2; }
[[ "$MAX" =~ ^[0-9]+$ ]] || { echo "[ERR] --max must be integer" >&2; exit 2; }
[[ "$MAPQ" =~ ^[0-9]+$ ]] || { echo "[ERR] -q must be integer" >&2; exit 2; }
[[ "$BASEQ" =~ ^[0-9]+$ ]] || { echo "[ERR] -Q must be integer" >&2; exit 2; }
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "[ERR] --threads must be integer" >&2; exit 2; }
(( MAX >= MIN )) || { echo "[ERR] --max (${MAX}) must be >= --min (${MIN})" >&2; exit 2; }

mkdir -p "$OUTDIR"

# -------------------- main --------------------

# Read list, skipping blanks and comments; strip CR if DOS-formatted.
while IFS= read -r bam; do
  bam="${bam%$'\r'}"
  [[ -z "$bam" || "$bam" =~ ^# ]] && continue

  if [[ ! -s "$bam" ]]; then
    echo "[WARN] Skipping missing BAM: $bam"
    continue
  fi

  sample="$(basename "${bam%.bam}")"
  mask="${OUTDIR%/}/${sample}.mask.bed"

  echo "[INFO] Processing ${sample}  (min=${MIN}, max=${MAX}, q=${MAPQ}, Q=${BASEQ}, threads=${THREADS})"

  # Ensure BAM index exists
  if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
    echo "[INFO] Indexing BAM: $bam"
    samtools index -@ "$THREADS" "$bam"
  fi

  # Compute depth (including zero-coverage sites), flag outliers, convert to BED, sort, merge.
  # samtools depth columns: contig  pos(1-based)  depth
  samtools depth -a -q "$MAPQ" -Q "$BASEQ" -@ "$THREADS" "$bam" \
    | awk -v L="$MIN" -v H="$MAX" '($3<L || $3>H){ start=$2-1; if(start<0) start=0; print $1"\t"start"\t"$2 }' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > "$mask"

  # Quick feedback
  nlines=$(wc -l < "$mask" || echo 0)
  echo "[OK] Wrote ${mask}  (merged intervals: ${nlines})"

done < "$LIST"

echo "[DONE] Coverage masks in ${OUTDIR}"

