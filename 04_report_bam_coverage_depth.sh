#!/usr/bin/env bash
set -euo pipefail

# report_bam_coverage_depth.sh
#
# Compute length-weighted mean coverage and depth histograms for BAMs.
#
# Usage:
#   bash report_bam_coverage_depth.sh sample1.bam sample2.bam
#   bash report_bam_coverage_depth.sh -f bam_list.txt
#   bash report_bam_coverage_depth.sh -f bam_list.txt --contigs contigs.txt
#
# Options:
#   -f FILE        File with BAM paths (one per line)
#   -o FILE        Output TSV for coverage (default: ./bam_coverage.tsv)
#   --contigs FILE Plain-text contig list (one contig name per line)
#   --regex REGEX  Filter contigs by regex (ignored if --contigs is used)
#   --histdir DIR  Output directory for depth histograms (default: ./depth_hists)
#   --label-col S  Header name for sample column (default: sample)
#
# Output:
#   - Coverage TSV:   <label-col>  bam  mean_coverage
#   - Histograms:     depth_hists/<sample>.depth.hist  (two columns: depth count)
#
# Notes:
#   - Assumes BAMs already filtered to primary alignments; supplementary is excluded (--ff SUPPLEMENTARY).
#   - If --contigs is provided, coverage is computed by filtering samtools output rows
#     to those contigs, and depth is computed by generating a BED from BAM header lengths.

LIST=""
OUT="$(pwd)/bam_coverage.tsv"
CONTIGS_FILE=""
REGEX=""
HISTDIR="./depth_hists"
LABEL_COL="sample"

usage() { sed -n '2,120p' "$0" | sed 's/^# \{0,1\}//'; exit 0; }

ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -f) LIST="$2"; shift 2 ;;
    -o) OUT="$2"; shift 2 ;;
    --contigs) CONTIGS_FILE="$2"; shift 2 ;;
    --regex) REGEX="$2"; shift 2 ;;
    --histdir) HISTDIR="$2"; shift 2 ;;
    --label-col) LABEL_COL="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) ARGS+=("$1"); shift ;;
  esac
done

INPUTS=()
if [[ -n "$LIST" ]]; then
  while IFS= read -r line; do
    [[ -z "$line" ]] && continue
    INPUTS+=("$line")
  done < "$LIST"
fi
INPUTS+=("${ARGS[@]}")
[[ ${#INPUTS[@]} -gt 0 ]] || { echo "[ERR] No BAMs provided."; exit 2; }

[[ -f "$OUT" ]] || echo -e "${LABEL_COL}\tbam\tmean_coverage" > "$OUT"
mkdir -p "$HISTDIR"

# Build a regex from a contig-name list: ^(c1|c2|c3)$
build_regex_from_list() {
  local list="$1"
  awk 'NF{gsub(/[][^$.|?*+(){}\\]/,"\\&"); a[++n]=$0} END{
        if(n){printf("^("); for(i=1;i<=n;i++){printf("%s%s", (i>1?"|":""), a[i])} printf(")$\n")}
      }' "$list"
}

# Create a BED (chrom 0-based-start 1-based-end) covering full contigs in list, using BAM header SQ LN
make_full_contig_bed_from_bam() {
  local bam="$1"
  local contigs_list="$2"
  local outbed="$3"
  # map contig -> 1 in an assoc array from list
  awk 'NF{a[$0]=1} END{for(k in a) print k}' "$contigs_list" | sort > "${outbed}.names"
  # extract SN and LN from BAM header
  samtools view -H "$bam" \
    | awk -F'\t' '$1=="@SQ"{sn="";ln="";
        for(i=1;i<=NF;i++){ if($i~/^SN:/) sn=substr($i,4); else if($i~/^LN:/) ln=substr($i,4) }
        if(sn!="" && ln!=""){ print sn"\t"ln }
      }' | sort > "${outbed}.all"
  # join to keep only requested contigs and write BED
  join -t $'\t' -1 1 -2 1 "${outbed}.names" "${outbed}.all" \
    | awk -F'\t' '{printf "%s\t0\t%d\n",$1,$2}' > "$outbed"
  rm -f "${outbed}.names" "${outbed}.all"
}

compute_cov() {
  local bam="$1" contigs_file="$2" regex="$3"
  if [[ -n "$contigs_file" ]]; then
    local rx; rx="$(build_regex_from_list "$contigs_file")"
    samtools coverage --ff SUPPLEMENTARY "$bam" \
    | awk -v rx="$rx" 'NR==1{next} $1 ~ rx {len=$3-$2+1; s+=$7*len; t+=len} END{if(t) printf "%.4f\n", s/t; else print "0.0000"}'
  elif [[ -n "$regex" ]]; then
    samtools coverage --ff SUPPLEMENTARY "$bam" \
    | awk -v rx="$regex" 'NR==1{next} $1 ~ rx {len=$3-$2+1; s+=$7*len; t+=len} END{if(t) printf "%.4f\n", s/t; else print "0.0000"}'
  else
    samtools coverage --ff SUPPLEMENTARY "$bam" \
    | awk 'NR==1{next} {len=$3-$2+1; s+=$7*len; t+=len} END{if(t) printf "%.4f\n", s/t; else print "0.0000"}'
  fi
}

compute_hist() {
  local bam="$1" contigs_file="$2" sample="$3" outdir="$4"
  local outfile="${outdir}/${sample}.depth.hist"
  if [[ -n "$contigs_file" ]]; then
    local bed="${outfile}.tmp.bed"
    make_full_contig_bed_from_bam "$bam" "$contigs_file" "$bed"
    # samtools depth expects BED (0-based start, 1-based end)
    samtools depth -a -b "$bed" "$bam" \
      | awk '{d[$3]++} END{for (k in d) print k,d[k]}' | sort -n > "$outfile"
    rm -f "$bed"
  else
    samtools depth -a "$bam" \
      | awk '{d[$3]++} END{for (k in d) print k,d[k]}' | sort -n > "$outfile"
  fi
  echo "[HIST] Wrote $outfile"
}

for bam in "${INPUTS[@]}"; do
  [[ -s "$bam" ]] || { echo "[WARN] Skipping missing BAM: $bam"; continue; }
  sample="$(basename "${bam%.bam}")"
  echo "[INFO] Processing $sample"

  mean_cov="$(compute_cov "$bam" "$CONTIGS_FILE" "$REGEX")"

  tmp="${OUT}.tmp.$$"
  awk -v s="$sample" -F'\t' 'NR==1{print;next} $1!=s' "$OUT" > "$tmp" || true
  echo -e "${sample}\t${bam}\t${mean_cov}" >> "$tmp"
  mv "$tmp" "$OUT"

  compute_hist "$bam" "$CONTIGS_FILE" "$sample" "$HISTDIR"
  echo "[STAT] ${sample}: mean coverage = ${mean_cov}×"
done

echo "[OK] Coverage summary: $OUT"
echo "[OK] Depth histograms in: $HISTDIR"

