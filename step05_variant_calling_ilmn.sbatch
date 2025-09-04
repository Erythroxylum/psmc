#!/bin/bash
#
#SBATCH -p sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH -t 06:00:00
#SBATCH -o log.%j.out
#SBATCH -e log.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<EMAIL>
# =========================
# EDIT SCRIPTDIR PATH PER SYSTEM
# =========================
# Hard-code SCRIPTDIR here:
SCRIPTDIR="$HOME/scripts/psmc/05_variant_calling"
#
#
# Usage examples:
#   sbatch --job-name=vcf-ILL step05_variant_calling_ilmn.sbatch -r ref.fa -f illumina_bams.txt
#   sbatch --job-name=vcf-ILL step05_variant_calling_ilmn.sbatch -r ref.fa sample1.illumina.singleton.bam sample2.illumina.singleton.bam
#
# This wraps 05_variant_calling_ilmn.sh


# Require at least one argument (forwarded to the script)
[[ $# -ge 1 ]] || { echo "Usage: $0 <args to 05_variant_calling_ilmn.sh>"; exit 2; }

# Path to script
SCRIPT="$SCRIPTDIR/05_variant_calling_ilmn.sh"

# Activate conda env
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate psmc

echo "[INFO] Host: $(hostname)"
echo "[INFO] Date: $(date)"

# Run
bash "$SCRIPT" "$@"

echo "[OK] Finished Illumina variant calling."

