#!/usr/bin/env bash
#SBATCH -p test
#SBATCH -J psmcBootMC
#SBATCH -o logs/bootmc_%j.out
#SBATCH -e logs/bootmc_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=2G             # memory for all concurrent procs
#SBATCH -t 03:00:00           # walltime
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<EMAIL>

set -euo pipefail

# -----------------------------
# Usage:
# sbatch --cpus-per-task=12 step07_psmc_boot_multicore.sbatch \
#   <input.psmcfa> <outdir> "<pattern>" <t> <r> <Niter> <NBOOT> [NPROCS]
#
# Example:
# sbatch --cpus-per-task=10 step07_psmc_boot_multicore.sbatch \
#   input.psmcfa psmc_pfurb/02_standard "4+25*2+4+6" 10 5 25 100 10
# -----------------------------

if [ "$#" -lt 7 ]; then
  echo "Usage: $0 <input.psmcfa> <outdir> \"<pattern>\" <t> <r> <Niter> <NBOOT> [NPROCS]" >&2
  exit 1
fi

IN_PSMCFA=$1
OUTDIR=$2
PATT=$3
TVAL=$4
RVAL=$5
NITER=$6
NBOOT=$7
NPROCS=${8:-${SLURM_CPUS_PER_TASK:-1}}

mkdir -p "${OUTDIR}/boot" logs

# Make/locate split.psmcfa (once)
SPLIT="${OUTDIR}/split.psmcfa"
if [ ! -s "${SPLIT}" ]; then
  if [ -s "${IN_PSMCFA}" ]; then
    ln -sf "$(readlink -f "${IN_PSMCFA}")" "${OUTDIR}/input.psmcfa"
    ( cd "${OUTDIR}" && splitfa input.psmcfa > split.psmcfa )
  else
    echo "ERROR: input PSMCFA not found: ${IN_PSMCFA}" >&2
    exit 1
  fi
fi

# Run main model once (skip if exists)
MAIN_OUT="${OUTDIR}/out_main.psmc"
if [ ! -s "${MAIN_OUT}" ]; then
  echo "Running main PSMC model..."
  psmc -p "${PATT}" -t "${TVAL}" -r "${RVAL}" -N "${NITER}" \
       -o "${MAIN_OUT}" "${SPLIT/\/split.psmcfa//input.psmcfa}"
fi

# Environment hygiene: keep each PSMC single-threaded
export OMP_NUM_THREADS=1

echo "Launching ${NBOOT} bootstraps with up to ${NPROCS} concurrent processes..."
# Prefer GNU parallel if available
if command -v parallel >/dev/null 2>&1; then
  parallel --jobs "${NPROCS}" \
    'psmc -b -p "'"${PATT}"'" -t '"${TVAL}"' -r '"${RVAL}"' -N '"${NITER}"' \
     -o "'"${OUTDIR}"'/boot/boot_{1}.psmc" "'"${SPLIT}"'"' ::: $(seq 1 "${NBOOT}")
else
  # Pure bash throttle (requires bash >=4.3 for wait -n)
  running=0
  for i in $(seq 1 "${NBOOT}"); do
    psmc -b -p "${PATT}" -t "${TVAL}" -r "${RVAL}" -N "${NITER}" \
         -o "${OUTDIR}/boot/boot_${i}.psmc" "${SPLIT}" &
    (( running++ ))
    if (( running >= NPROCS )); then
      wait -n
      (( running-- ))
    fi
  done
  wait
fi

echo "All bootstraps finished for ${OUTDIR}"

