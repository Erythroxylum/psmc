#!/usr/bin/env bash
set -euo pipefail
N=$(wc -l < fq.list)

# Common scaling for P. furbishiae
MU=7.0e-9
GEN=3
Q=20
GMIN=10000

# sbatch --cpus-per-task=12 step07_psmc_boot_multicore.sbatch \
#   <input.psmcfa> <outdir> "<pattern>" <t> <r> <Niter> <NBOOT> [NPROCS]

# Conservative 
sbatch --cpus-per-task=10 psmc_boot_multicore.sbatch \
 psmc_out_pfurb/01_conservative/Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton.iupac.min500000/Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton.iupac.min500000.split.psmcfa \
 psmc_out_pfurb/01_conservative "18+16*2+4+4" 10 5 25 100 10 

# Standard 
sbatch --cpus-per-task=10 psmc_boot_multicore.sbatch \
 psmc_out_pfurb/02_standard/Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton.iupac.min500000.split.psmcfa \
 psmc_out_pfurb/02_standard "4+25*2+4+6" 10 5 25 100 10 

# Recent-focused 
sbatch --cpus-per-task=20 psmc_boot_multicore.sbatch \
 psmc_out_pfurb/Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton.iupac.min500000.split.psmcfa \
 psmc_out_pfurb/03_recent_focus "10*1+14*2+4*4+10" 10 5 30 100 20 

# Aggressive recent 
sbatch --cpus-per-task=20 psmc_boot_multicore.sbatch \
 psmc_out_pfurb/Pedicularis_furbishiae_Pfurb.ill.sorted.markdup.singleton.iupac.min500000.split.psmcfa \
 psmc_out_pfurb/04_aggressive_recent "15*1+10*2+4*4+13" 8 5 30 100 20

