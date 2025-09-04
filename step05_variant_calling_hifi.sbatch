#!/bin/bash
#
#SBATCH -p test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o log.%j.out
#SBATCH -e log.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<EMAIL>
#
# =========================
# EDIT SCRIPTDIR PATH PER SYSTEM
# =========================
# Hard-code SCRIPTDIR here:
SCRIPTDIR="$HOME/scripts/psmc/05_variant_calling"
#
#
# Usage examples:
#   sbatch --job-name=vcf-HIFI step05_variant_calling_hifi.sbatch -r ref.fa -f hifi_bams.txt
#   sbatch --job-name=vcf-HIFI step05_variant_calling_hifi.sbatch -r ref.fa sample1.hifi.nuclear.bam sample2.hifi.nuclear.bam
#
# This wraps 05_variant_calling_hifi.sh

[[ $# -ge 1 ]] || { echo "Usage: $0 <args to 05_variant_calling_hifi.sh>"; exit 2; }

SCRIPT="$SCRIPTDIR/05_variant_calling_hifi.sh"

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate psmc

echo "[INFO] Host: $(hostname)"
echo "[INFO] Date: $(date)"

bash "$SCRIPT" "$@"

echo "[OK] Finished HiFi variant calling."

