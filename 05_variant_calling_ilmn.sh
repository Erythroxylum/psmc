#!/usr/bin/env bash
set -euo pipefail

# 05_variant_calling_ilmn.sh
#
# Variant calling for Illumina BAMs using bcftools mpileup + call.
#
# Usage:
#   bash 05_variant_calling_ilmn.sh -r REF.fa -f bam_list.txt
#   bash 05_variant_calling_ilmn.sh -r REF.fa /path/S1.illumina.singleton.bam /path/S2.illumina.singleton.bam
#
# Required:
#   -r, --ref REF.fa      Reference FASTA (indexed: REF.fa.fai)
#
# Optional:
#   -f, --bam-list FILE   File with BAM paths (one per line)
#   -o, --outdir DIR      Output directory for VCFs (default: ./vcf_ilmn)
#   -q INT                Min mapping quality for mpileup (default: 20)
#   -Q INT                Min base quality for mpileup (default: 20)
#   -C INT                Coefficient for downgrading MQ near indels (default: 50)
#   --contigs FILE        Optional list of contig names to restrict calling (one per line)
#   --set-sample          Rename VCF sample to the BAM basename (uses bcftools reheader)
#
# Output (per BAM):
#   <outdir>/<sample>.ilmn.vcf.gz and .tbi
#
# Notes:
#   * Assumes BAMs are already primary-only (no unmapped/secondary/supplementary/dup/QC-fail).
#   * Depth outlier masking is best applied later at the consensus step (bcftools consensus -m mask.bed).
#   * If you supply --contigs, make sure names match the BAM/REF header exactly.

REF=""
BAMLIST=""
OUTDIR="./vcf_ilmn"
MAPQ=20
BASEQ=20
COEFF=50
CONTIGS_FILE=""
SET_SAMPLE=0

# --- parse args ---
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--ref) REF="${2:?}"; shift 2 ;;
    -f|--bam-list) BAMLIST="${2:?}"; shift 2 ;;
    -o|--outdir) OUTDIR="${2:?}"; shift 2 ;;
    -q) MAPQ="${2:?}"; shift 2 ;;
    -Q) BASEQ="${2:?}"; shift 2 ;;
    -C) COEFF="${2:?}"; shift 2 ;;
    --contigs) CONTIGS_FILE="${2:?}"; shift 2 ;;
    --set-sample) SET_SAMPLE=1; shift 1 ;;
    -h|--help) sed -n '2,160p' "$0"; exit 0 ;;
    *) ARGS+=("$1"); shift ;;
  esac
done

[[ -n "$REF" && -s "$REF" ]] || { echo "[ERR] Provide -r REF.fa"; exit 2; }
[[ -s "${REF}.fai" ]] || { echo "[ERR] Missing ${REF}.fai (run: samtools faidx ${REF})"; exit 2; }

# tool checks we rely on
command -v bcftools >/dev/null 2>&1 || { echo "[ERR] bcftools not found"; exit 3; }
command -v tabix     >/dev/null 2>&1 || { echo "[ERR] tabix (htslib) not found"; exit 3; }

INPUTS=()
if [[ -n "$BAMLIST" ]]; then
  while IFS= read -r line; do [[ -n "$line" ]] && INPUTS+=("$line"); done < "$BAMLIST"
fi
INPUTS+=("${ARGS[@]}")
[[ ${#INPUTS[@]} -gt 0 ]] || { echo "[ERR] No BAMs provided."; exit 2; }

mkdir -p "$OUTDIR"

# Build -r args if contigs file supplied (mpileup accepts multiple -r)
RARGS=()
if [[ -n "$CONTIGS_FILE" ]]; then
  while IFS= read -r c; do [[ -n "$c" ]] && RARGS+=(-r "$c"); done < "$CONTIGS_FILE"
fi

for BAM in "${INPUTS[@]}"; do
  [[ -s "$BAM" ]] || { echo "[WARN] Skipping missing BAM: $BAM"; continue; }
  SAMPLE="$(basename "${BAM%.bam}")"
  OUTVCF="${OUTDIR}/${SAMPLE}.ilmn.vcf.gz"

  echo "[INFO] Calling (Illumina): $SAMPLE"
  # mpileup -> call
  bcftools mpileup -f "$REF" -q "$MAPQ" -Q "$BASEQ" -C "$COEFF" "${RARGS[@]}" "$BAM" \
  | bcftools call -mv \
  | bgzip -c > "$OUTVCF"

  tabix -p vcf -f "$OUTVCF"

  # Optional: set VCF sample name to BAM basename
  if [[ "$SET_SAMPLE" -eq 1 ]]; then
    smap="$(mktemp)"
    echo "$SAMPLE" > "$smap"
    out2="${OUTDIR}/${SAMPLE}.ilmn.sample.vcf.gz"
    bcftools reheader -s "$smap" -o "$out2" "$OUTVCF"
    tabix -f -p vcf "$out2"
    mv -f "$out2" "$OUTVCF"
    mv -f "$out2.tbi" "$OUTVCF.tbi"  # rename index properly
    rm -f "$smap"
    echo "[INFO] Renamed VCF sample to: $SAMPLE"
  fi

  echo "[OK] $OUTVCF"
done

echo "[DONE] Illumina variant calling complete → $OUTDIR"

