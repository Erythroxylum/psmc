#!/usr/bin/env bash
set -euo pipefail

# 05_variant_calling_hifi.sh
#
# Variant calling for PacBio HiFi BAMs using Longshot (v1.0.0).
#
# Usage:
#   # BAMs directly
#   bash 05_variant_calling_hifi.sh -r REF.fa /path/S1.hifi.bam /path/S2.hifi.bam
#
#   # BAM list file
#   bash 05_variant_calling_hifi.sh -r REF.fa -f hifi_bams.txt
#
# Optional:
#   --contigs FILE   Plain-text list of contig names (one per line) to restrict calling.
#                    If provided, the script calls Longshot per contig with `-r <contig>` and concatenates results.
#   -o, --outdir D   Output directory (default: ./vcf_hifi)
#   -s, --set-sample Use BAM basename (without .bam) as VCF sample ID via Longshot `-s`. (off by default)
#
# Notes:
#   * Requires: longshot (>=1.0.0), bcftools, bgzip/tabix (htslib), samtools.
#   * This script does NOT pass any non-existent flags to longshot (no threads flag).
#   * Safe defaults: we add Longshot `-F` (force overwrite) so re-runs clobber temp files.
#   * Depth/region masking should be applied later during consensus (e.g., bcftools consensus -m mask.bed).
#
# Output (per BAM):
#   <outdir>/<sample>.hifi.vcf.gz  and  <outdir>/<sample>.hifi.vcf.gz.tbi
#
# Exit codes:
#   2: missing required args; 3: missing tools; 4: missing reference index

REF=""
BAMLIST=""
OUTDIR="./vcf_hifi"
CONTIGS_FILE=""
SET_SAMPLE_ID=0

# --- parse args ---
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--ref) REF="${2:?}"; shift 2 ;;
    -f|--bam-list) BAMLIST="${2:?}"; shift 2 ;;
    -o|--outdir) OUTDIR="${2:?}"; shift 2 ;;
    --contigs) CONTIGS_FILE="${2:?}"; shift 2 ;;
    -s|--set-sample) SET_SAMPLE_ID=1; shift 1 ;;
    -h|--help) sed -n '2,200p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) ARGS+=("$1"); shift ;;
  esac
done

# --- required args ---
[[ -n "$REF" && -s "$REF" ]] || { echo "[ERR] Provide -r REF.fa"; exit 2; }
[[ -s "${REF}.fai" ]] || { echo "[ERR] Missing ${REF}.fai (run: samtools faidx ${REF})"; exit 4; }

# --- tool checks ---
need() { command -v "$1" >/dev/null 2>&1 || { echo "[ERR] Missing tool: $1"; exit 3; }; }
need longshot
need bcftools
need bgzip
need tabix
need samtools

# --- build BAM input list ---
INPUTS=()
if [[ -n "$BAMLIST" ]]; then
  while IFS= read -r line; do [[ -n "$line" ]] && INPUTS+=("$line"); done < "$BAMLIST"
fi
INPUTS+=("${ARGS[@]}")
[[ ${#INPUTS[@]} -gt 0 ]] || { echo "[ERR] No BAMs provided."; exit 2; }

# --- prep ---
mkdir -p "$OUTDIR"

# --- core: call with longshot, optionally per contig, then bgzip+tabix ---
call_longshot() {
  local bam="$1" ref="$2" outvcf="$3" set_sid="$4" contigs_file="${5:-}"
  local sample_id
  sample_id="$(basename "${bam%.bam}")"

  if [[ -n "$contigs_file" ]]; then
    # Per-contig calling, then concatenate
    local tmpdir; tmpdir="$(mktemp -d)"
    trap 'rm -rf "$tmpdir"' RETURN

    local parts=()
    while IFS= read -r chr; do
      [[ -z "$chr" ]] && continue
      local part="${tmpdir}/${sample_id}.${chr}.vcf"
      if [[ "$set_sid" -eq 1 ]]; then
        longshot -b "$bam" -f "$ref" -r "$chr" -o "$part" -F -s "$sample_id"
      else
        longshot -b "$bam" -f "$ref" -r "$chr" -o "$part" -F
      fi
      parts+=("$part")
    done < "$contigs_file"

    # bgzip each part and concat
    local zparts=()
    for p in "${parts[@]}"; do
      bgzip -f "$p"
      zparts+=("${p}.gz")
    done

    # Safe concat; bcftools handles headers
    bcftools concat -Oz -o "$outvcf" "${zparts[@]}"

  else
    # Whole BAM in one shot
    local tmpvcf="${outvcf%.gz}.tmp.vcf"
    if [[ "$set_sid" -eq 1 ]]; then
      longshot -b "$bam" -f "$ref" -o "$tmpvcf" -F -s "$sample_id"
    else
      longshot -b "$bam" -f "$ref" -o "$tmpvcf" -F
    fi
    bgzip -f "$tmpvcf"
    mv "${tmpvcf}.gz" "$outvcf"
  fi

  tabix -f -p vcf "$outvcf"
}

# --- iterate over BAMs ---
for BAM in "${INPUTS[@]}"; do
  [[ -s "$BAM" ]] || { echo "[WARN] Skipping missing BAM: $BAM"; continue; }
  SAMPLE="$(basename "${BAM%.bam}")"
  OUTVCF="${OUTDIR}/${SAMPLE}.hifi.vcf.gz"

  echo "[INFO] Longshot calling: $SAMPLE"
  call_longshot "$BAM" "$REF" "$OUTVCF" "$SET_SAMPLE_ID" "$CONTIGS_FILE"
  echo "[OK] $OUTVCF"
done

echo "[DONE] HiFi variant calling complete → $OUTDIR"

