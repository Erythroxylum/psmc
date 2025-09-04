#!/bin/bash
#SBATCH -p sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
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
# Hard-code SCRIPTDIR here.
SCRIPTDIR="$HOME/scripts/psmc/03_map"
#
# Usage:
#   sbatch --job-name=map-HIFI-<TAG> step03_map_hifi.sbatch <TAG> <REF> <HiFi.clean.gz> <PLATFORM>



[[ $# -eq 4 ]] || { echo "Usage: $0 <TAG> <REF> <HiFi.clean.gz> <PLATFORM>"; exit 2; }
TAG="$1"; REF="$2"; READS="$3"; PLATFORM="$4"

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate psmc

bash "$SCRIPTDIR/03_map_hifi.sh" "$TAG" "$REF" "$READS" "$PLATFORM"

