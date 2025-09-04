#!/bin/bash
#
# Usage:
#   sbatch --job-name=hifi-ERRxxxxxxx step02_hifi.sbatch ERRxxxxxxx
#
# Example:
#   sbatch --job-name=hifi-ERR12370312 step02_hifi.sbatch ERR12370312
#
# Input:
#   BASEDIR/reads/[RUN].fastq(.gz)
# Output:
#   BASEDIR/clean_reads/[RUN].min1000.q25.p90.hifi.fastq.gz

#SBATCH -p sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o log.%j.out
#SBATCH -e log.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<EMAIL>

# --- CONFIG ---
if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch --job-name=hifi-RUN $0 RUN" >&2
    exit 1
fi
RUN="$1"

# --- ENV ---
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate psmc

# --- RUN ---
bash 02_hifi_clean.sh "$RUN"

