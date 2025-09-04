#!/usr/bin/env bash
#SBATCH -p sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o log.%j.out
#SBATCH -e log.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<EMAIL>
#
# =========================
# EDIT SCRIPTDIR PATH PER SYSTEM
# =========================
# Hard-code SCRIPTDIR here:
SCRIPTDIR="$HOME/scripts/psmc/03_map"
#
# -------------------------------------------------------------------------
# Usage examples:
#
# 1. Just sort & index BAM (no filtering):
#    sbatch step03b_filter_bam.sbatch sample.bam
#
# 2. Filter Illumina BAMs to primary, mapped, non-dup, non-QC-fail:
#    sbatch step03b_filter_bam.sbatch --exclude-flags 3844 sample1.bam sample2.bam
#
# 3. Subset to specific contigs (any platform):
#    sbatch step03b_filter_bam.sbatch --contigs chr.list sample.bam
#
# 4. Use a list of BAMs:
#    sbatch step03b_filter_bam.sbatch -f bam_list.txt --exclude-flags 2308
#
# Notes:
# - Outputs are written next to inputs as <in>.filtered.bam (+ .bai).
# - This script is platform-agnostic; works for Illumina or HiFi.
# - Filtering is explicit-only: no contigs or flags are applied unless given.
# -------------------------------------------------------------------------

set -euo pipefail

# Resolve worker script path
WORKER="${SCRIPTDIR}/03b_filter_bam.sh"

if [[ ! -x "$WORKER" ]]; then
  echo "[ERR] Cannot find executable 03_filter_bam.sh at: $WORKER" >&2
  exit 2
fi

bash "$WORKER" "$@"
