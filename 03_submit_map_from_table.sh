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
# Run full pipeline from a samples table.
#
# Usage:
#   bash 03_submit_map_from_table.sh /path/to/samples_psmc.txt [TARGET_COV] [--reset-stats]
#
# No downsampling:
# bash ~/scripts/psmc/map/03_submit_map_from_table.sh ~/scripts/psmc/samples_psmc.txt
#
# Downsample all mappings to 20× with rasusa:
# bash ~/scripts/psmc/map/03_submit_map_from_table.sh ~/scripts/psmc/samples_psmc.txt 10
#
# Fresh header + 3× test run:
# bash ~/scripts/psmc/map/03_submit_map_from_table.sh ~/scripts/psmc/samples_psmc.txt 3 --reset-stats
#
#
# Notes:
# - TARGET_COV is optional. If provided, it is exported to mapping jobs to enable rasusa downsampling.
# - Use --reset-stats to recreate map/map_stats.txt header before submissions.
#
# Table columns (tab-delimited), expected headers incl. (order per your file):
# Assembly Accession, Assembly Name, BioProject, BioSample, SRA, seq_platform, sci_name, genome_fasta,
# read1_fastq, read2_fastq, clean_read1_fastq, clean_read2_fastq, ... (rest ignored)


TABLE="${1:?Provide samples_psmc.txt}"
TARGET_COV_ARG="${2:-}"
RESET_FLAG="${3:-}"

SCRIPTDIR="$HOME/scripts/psmc"
MAP_DIR="${BASEDIR}/map"
STATS_FILE="${MAP_DIR}/map_stats.txt"

mkdir -p "${MAP_DIR}"

# Optionally (re)initialize map_stats.txt header
if [[ "${RESET_FLAG:-}" == "--reset-stats" ]] || [[ ! -f "${STATS_FILE}" ]]; then
  echo -e "sample\tplatform\tbam\treads_total\treads_mapped\tpct_mapped\tmean_coverage\ttarget_cov\tachieved_cov" > "${STATS_FILE}"
  echo "[INIT] Wrote header to ${STATS_FILE}"
fi

# Helper: submit sbatch and return jobid
submit() { sbatch --parsable "$@"; }

# Optional downsampling export
EXPORT_OPT=()
if [[ -n "${TARGET_COV_ARG}" ]]; then
  # ensure numeric-ish
  if [[ "${TARGET_COV_ARG}" =~ ^[0-9]+$|^[0-9]+(\.[0-9]+)?$ ]]; then
    EXPORT_OPT=(--export=ALL,TARGET_COV="${TARGET_COV_ARG}")
    echo "[INFO] Will export TARGET_COV=${TARGET_COV_ARG} to mapping jobs"
  else
    echo "[WARN] TARGET_COV '${TARGET_COV_ARG}' not numeric; ignoring"
  fi
fi

# Print a header for job tracking
printf "%-35s %-12s %-12s %-14s %-14s %-14s\n" "TAG" "PLATFORM" "SRA" "JOB_FETCH" "JOB_CLEAN" "JOB_MAP"
printf "%-35s %-12s %-12s %-14s %-14s %-14s\n" "-----------------------------------" "--------" "------------" "------------" "------------" "------------"

# Parse table skipping header, extracting needed fields
# Fields (1-based): 5=SRA, 6=seq_platform, 7=sci_name, 8=genome_fasta,
# 9=read1_fastq, 10=read2_fastq, 11=clean_read1_fastq, 12=clean_read2_fastq
tail -n +2 "${TABLE}" | awk -F'\t' 'NF>0' | while IFS=$'\t' read -r col1 col2 col3 col4 SRA PLATFORM SCI REF R1 R2 C1 C2 rest; do
  # Skip blank lines
  [[ -z "${SCI// }" && -z "${SRA// }" ]] && continue

  # Normalize platform string
  platform_lc="$(echo "${PLATFORM:-}" | tr '[:upper:]' '[:lower:]')"
  if [[ "${platform_lc}" != "illumina" && "${platform_lc}" != "hifi" && "${platform_lc}" != "pacbio" ]]; then
    echo "[WARN] Unknown platform '${PLATFORM}' for ${SCI} (${SRA}); skipping."
    continue
  fi
  [[ "${platform_lc}" == "pacbio" ]] && platform_lc="hifi"

  # Build a tag: sci_name + SRA (or sci_name only if SRA missing)
  TAG="${SCI}"
  [[ -n "${SRA:-}" && "${SRA}" != "NA" ]] && TAG="${SCI}_${SRA}"

  # Trim whitespace
  TAG="${TAG// /_}"

  # Decide inputs and which scripts to run
  job_fetch="NA"
  job_clean="NA"
  job_map="NA"

  # Reference must exist
  if [[ ! -s "${REF}" ]]; then
    echo "[WARN] Missing REF for ${TAG}: ${REF}; skipping."
    continue
  fi

  if [[ "${platform_lc}" == "illumina" ]]; then
    # Prefer clean reads if present
    if [[ -s "${C1:-}" && -s "${C2:-}" ]]; then
      R1C="${C1}"; R2C="${C2}"
    else
      # If raw reads missing, try to fetch via SRA if given
      if [[ ! -s "${R1:-}" || ! -s "${R2:-}" ]]; then
        if [[ -n "${SRA:-}" && "${SRA}" != "NA" ]]; then
          job_fetch=$(submit --job-name="fetch-${SRA}" "${SCRIPTDIR}/step01_fetch.sbatch" "${SRA}")
        else
          echo "[WARN] No raw reads and no SRA for ${TAG}; skipping."
          continue
        fi
      fi
      # Clean (fastp) depends on fetch if fetch was submitted
      if [[ -n "${SRA:-}" && "${SRA}" != "NA" ]]; then
        dep=()
        [[ "${job_fetch}" != "NA" ]] && dep=(--dependency=afterok:"${job_fetch}")
        job_clean=$(submit "${dep[@]}" --job-name="fastp-${SRA}" "${SCRIPTDIR}/step02_fastp.sbatch" "${SRA}")
        # Paths produced by fastp.sh:
        R1C="/n/netscratch/davis_lab/Everyone/dwhite/psmc/clean_reads/${SRA}_R1.clean.fastq.gz"
        R2C="/n/netscratch/davis_lab/Everyone/dwhite/psmc/clean_reads/${SRA}_R2.clean.fastq.gz"
      else
        # No SRA; use given R1/R2 directly (assume already clean)
        R1C="${R1}"; R2C="${R2}"
      fi
    fi

    # Submit mapping; depend on clean if it ran
    dep=()
    [[ "${job_clean}" != "NA" ]] && dep=(--dependency=afterok:"${job_clean}")
    job_map=$(submit "${EXPORT_OPT[@]}" "${dep[@]}" \
      --job-name="map-ILL-${SRA:-${TAG}}" \
      "${SCRIPTDIR}/map/step03_map_illumina.sbatch" \
      "${TAG}" "${REF}" "${R1C}" "${R2C}" "illumina")

    printf "%-35s %-12s %-12s %-14s %-14s %-14s\n" "${TAG}" "illumina" "${SRA:-NA}" "${job_fetch}" "${job_clean}" "${job_map}"

  else
    # HiFi (single-end)
    # Prefer clean if present in C1 (we store HiFi clean in clean_read1_fastq)
    if [[ -s "${C1:-}" ]]; then
      HIFI_C="${C1}"
    else
      # If raw missing, fetch from SRA/ERR if available
      if [[ ! -s "${R1:-}" ]]; then
        if [[ -n "${SRA:-}" && "${SRA}" != "NA" ]]; then
          job_fetch=$(submit --job-name="fetch-${SRA}" "${SCRIPTDIR}/step01_fetch.sbatch" "${SRA}")
        else
          echo "[WARN] No HiFi reads and no SRA for ${TAG}; skipping."
          continue
        fi
      fi
      # Clean HiFi (filtlong) depends on fetch if it ran
      if [[ -n "${SRA:-}" && "${SRA}" != "NA" ]]; then
        dep=()
        [[ "${job_fetch}" != "NA" ]] && dep=(--dependency=afterok:"${job_fetch}")
        job_clean=$(submit "${dep[@]}" --job-name="hifi-${SRA}" "${SCRIPTDIR}/step02_hifi.sbatch" "${SRA}")
        HIFI_C="/n/netscratch/davis_lab/Everyone/dwhite/psmc/clean_reads/${SRA}.min1000.q25.p90.hifi.fastq.gz"
      else
        HIFI_C="${R1}"  # assume already clean if provided directly
      fi
    fi

    dep=()
    [[ "${job_clean}" != "NA" ]] && dep=(--dependency=afterok:"${job_clean}")
    job_map=$(submit "${EXPORT_OPT[@]}" "${dep[@]}" \
      --job-name="map-HIFI-${SRA:-${TAG}}" \
      "${SCRIPTDIR}/map/step03_map_hifi.sbatch" \
      "${TAG}" "${REF}" "${HIFI_C}" "hifi")

    printf "%-35s %-12s %-12s %-14s %-14s %-14s\n" "${TAG}" "hifi" "${SRA:-NA}" "${job_fetch}" "${job_clean}" "${job_map}"
  fi

done

