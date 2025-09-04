#!/bin/bash
#
# Usage:
#   sbatch --job-name=fastp-SRRxxxxxxx fastp-template.sbatch SRRxxxxxxx
#
# Example:
#   sbatch --job-name=fastp-SRR14835741 fastp-template.sbatch SRR14835741
#
# This will run fastp.sh on the paired reads in directory set in that script:
#   BASEDIR/reads/[SRR]_1.fastq
#   BASEDIR/reads/[SRR]_2.fastq
# and write cleaned reads + reports to:
#   BASEDIR/clean_reads/

#SBATCH -p sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o log.%j.out
#SBATCH -e log.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<EMAIL>

# --- CONFIG ---
if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch --job-name=fastp-SRRxxxxxxx $0 SRRxxxxxxx" >&2
    exit 1
fi
SRR="$1"

# --- ENV ---
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate psmc

# --- RUN ---
bash 02_fastp.sh "$SRR"
