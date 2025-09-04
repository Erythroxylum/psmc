#!/bin/bash
#
# Usage:
#   sbatch --job-name=fetch-SRRxxxxxxx step01_fetch.sbatch SRRxxxxxxx
#
# Example:
#   sbatch --job-name=fetch-SRR14835741 step01_fetch.sbatch SRR14835741
#
# This calls 01_fetch_reads_sra.sh, which downloads to outdir set in that script.
#  

#SBATCH -p sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o log.%j.out
#SBATCH -e log.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<EMAIL>

# --- args ---
if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch --job-name=fetch-RUN $0 RUN" >&2
  exit 1
fi
RUN="$1"

# --- env ---
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate psmc

# --- run ---
bash "$(dirname "$0")/01_fetch_reads_sra.sh" "$RUN"

