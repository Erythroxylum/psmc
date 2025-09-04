# PSMC Data Processing Pipeline

This project contains the scripts used to fetch raw reads, clean them, map to their genome assemblies, and prepare for PSMC analysis.

## Overview of Steps

#### Usage examples are provided in the header of each sbatch file. 

#### Users will also need to adjust slurm scripts (e.g., '-p') and output scratch directories for local servers in each shell and sbatch script.

## 0. Environment Setup
- **Script:** `00_make_psmc_env.sh`  
- **Purpose:** Create a `mamba` environment (`psmc`) with all required tools from the environment file `psmc_env.yml`. Environment is then activated.
  

## 1. Fetch Raw Reads
- **Script:** `01_fetch_reads_sra.sh` (worker)  
- **SBATCH Wrapper:** `step01_fetch.sbatch`  
- **Purpose:** Download SRA/ENA runs to project-specific directory specified in 01_fetch_reads_sra.sh. The script Download raw reads from SRA/ENA into a project-specific directory. The script (1) creates output and cache directories on scratch, (2) configures the SRA cache with vdb-config, (3) sets TMPDIR on scratch so temporary files don’t land in $HOME, (4) downloads runs with prefetch, converts to FASTQ with fasterq-dump, using --temp if supported.

## 2. Clean Reads
Illumina
- **Script:** `02_fastp.sh` (worker)
- **SBATCH Wrapper:** `step02_fastp.sbatch`
- **Purpose:** Adapter trim, quality filter, poly-G/X trim, deduplication.
- **Inputs:**  `reads/<RUN_ID>_1.fastq`, `reads/<RUN_ID>_2.fastq`
- **Outputs:** `clean_reads/*_R1.clean.fastq.gz`, `clean_reads/*_R2.clean.fastq.gz`.

PacBio HiFi Reads
- **Script:** `02_hifi_clean.sh` (worker)
- **SBATCH Wrapper:** `step02_hifi.sbatch`
- **Purpose:** Length/quality filter HiFi reads using `filtlong` (min length=1000 bp, min Q=25, keep top 90% by score).
- **Inputs:** `reads/<RUN_ID>.fastq` or ``.fastq.gz`
- **Outputs:** `clean_reads/<RUN_ID>.min1000.q25.p90.hifi.fastq.gz`.

## 3a. Mapping to References
  - Before mapping, download the reference assembly FASTA and remove organellar genomes or other unwanted scaffolds if necessary. 

  - Both Illumina and PacBio HiFi workflows will write mapping stats into `map/map_stats.txt`with the following tab-delimited columns:
`sample\tplatform\tbam\treads_total\treads_mapped\tpct_mapped\tmean_coverage\ttarget_cov\tachieved_cov`

- **Inputs:** Reference genome FASTA (per sample), cleaned reads.

#### Batch Submission from Table (uniform target coverage)
- **Script:** `03_submit_map_from_table.sh`
- **Purpose:** Read a tab-delimited sample table (`samples_psmc.txt`) and submit fetch, clean, and map jobs automatically.
- **Platforms supported:** `illumina` or (PacBio) `hifi` 
- **Target coverage:** A single `TARGET_COV` can be supplied and will be applied to all mapping jobs for psmc consistency. 20x-30x is ideal but this should be overshot if samples have contaminant reads.
- **Options:** `--reset-stats` resets the `map_stats.txt` header before submission.
- **Sample table columns (tab-delimited):**
`Assembly Accession, Assembly Name, BioProject, BioSample, SRA, seq_platform, sci_name, genome_fasta, read1_fastq, read2_fastq, clean_read1_fastq, clean_read2_fastq`

#### Independent Submissions (per-sample coverage)
  - Samples can also be submitted independently to specify a different TARGET_COV per sample. This may be necessary because contaminant reads that do not map reduce the achieved coverage. Per-sample scripts let you tune coverage individually.

Example for Illumina, 35x target:
`sbatch --export=ALL,TARGET_COV=35 --time=16:00:00 --cpus-per-task=8 \
  --job-name=map-DD \
  step03_map_illumina_bwa-mem2.sbatch \
  Dryas_drummondii_SRR5313982 \
  /path/to/ref.fa \
  /path/to/SRR5313982_R1.clean.fastq.gz \
  /path/to/SRR5313982_R2.clean.fastq.gz \
  illumina`
  
Example for HiFi, 35x target:
`sbatch --export=ALL,TARGET_COV=35 --time=16:00:00 --cpus-per-task=8 \
  --job-name=map-HIFI-Dryoct \
  step03_map_hifi.sbatch \
  Dryas_octopetala_ERR12370312 \
  /path/to/ref.fa \
  /path/to/ERR12370312.min1000.q25.p90.hifi.fastq.gz \
  hifi`

Illumina:
- **Script:** `03_map_illumina_bwa-mem2.sh`
- **SBATCH Wrapper:** `step03_map_illumina_bwa-mem2.sbatch`
- **Aligner:** `bwa-mem2 mem`
- **Optional:** Downsample with `rasusa` using `TARGET_COV`
- **Post-processing:** `samtools fixmate`, `samtools sort`, `samtools markdup`
- **Outputs:** `<SAMPLE>.ill.sorted.markdup.bam` + `bai`, Stats row in `map_stats.txt`.

HiFi:
- **Script:** `03_map_hifi.sh`
- **SBATCH Wrapper:** `step03_map_hifi.sbatch`
- **Aligner:** `minimap2 -ax map-hifi`, sort, index.
- **Optional:** Downsample with `rasusa` using `TARGET_COV`
- **Outputs:** `<SAMPLE>.hifi.sorted.bam` + `bai`, Stats row in `map_stats.txt`.


### 3b. Filter bams
  - Filtering is handled by a single, platform-agnostic script. By default, it does no filtering — it simply sorts and indexes BAMs. 
  - Filtering is only applied when explicitly requested via options.

Illumina
- **Script:** `03b_filter_bam.sh`
- **SBATCH Wrapper:** `step03b_filter_bam.sbatch`
- **Purpose:** 
  - Explicitly filter BAMs by flags (--exclude-flags, --include-flags) to drop or require specific SAM flags, or contigs (--contigs <file>) to keep only specified chromosomes.
  - If no filters are provided, the script just produces a coordinate-sorted, indexed copy of the input BAM. 
  - To keep only high-quality primary alignments; excludes unmapped, secondary, supplementary, QC-fail, and duplicate reads, use flag mask `-F 3844`. 
  - Script will print `samtools idxstats | head -10` for input and output bam files.
- **Inputs:** One or more mapped BAMs (`<SAMPLE>.ill.sorted.markdup.bam`) passed as arguments or as a list (`-f bam_list.txt`).
- **Outputs:** For each input `<SAMPLE>.bam`, writes `<SAMPLE>.singleton.bam` + `bai` in place.


## 4. Coverage Reporting and Depth Histograms

- **Script:** `04_report_bam_coverage_depth.sh`  
- **Inputs:** One or more **filtered BAMs** (produced in Step 3b). BAMs should already contain only primary, high-quality alignments.  
- **Purpose:**  
  - Compute length-weighted **mean coverage** per BAM using `samtools coverage`.  
  - Generate **depth histograms** (`depth vs. count of sites`) using `samtools depth`.  
  - Optionally restrict analysis to:  
    - A provided list of contigs (`--contigs file`)  
    - A regular expression filter (`--regex pattern`)  
- **Outputs:**  
  - A tab-delimited summary file `bam_coverage.tsv` with columns:  
    ```
    sample    bam    mean_coverage
    ```  
  - Depth histograms for each BAM in `depth_hists/`, e.g. `sample.depth.hist` (two columns: depth, count).  
- **Use:**  
  - This step is intended to be run **after filtering BAMs** to verify that desired coverage has been achieved.  
  - Examine the histograms to identify appropriate coverage cutoffs (`CLOW`/`CHIGH`, e.g. 10–60).  
  - These cutoffs and masks will later be used in **Step 6** to filter contigs by depth or exclude problematic regions.  


## 5. Variant Calling
- Both pipelines assume BAMs have already been filtered (Step 3b).  
- Masking of depth outliers or other BED-based filters is applied later, during consensus generation (Step 6). 

Illumina
- **Script:** `05_variant_calling_ilmn.sh`
- **SBATCH wrapper:** `step05_variant_calling_ilmn.sbatch`
- **Caller:** `bcftools mpileup` + `bcftools call`
- **Inputs:** Filtered Illumina BAMs and indexed reference FASTA (`REF.fa` + `REF.fa.fai`).
- **Options:** 
  - `-q`, `-Q`, `-C` for mapping/base quality thresholds and MQ downgrading (defaults: 20/20/50).  
  - `--contigs FILE` to restrict variant calling to a provided list of contigs.  
  - `--set-sample` to rename the VCF sample to match the BAM basename.
- **Outputs:** Per-BAM VCFs in `vcf_ilmn/` (default directory):  
  `<sample>.ilmn.vcf.gz` + index `.tbi`
- **Effort:** For Illumina, a small genome (~256 Mbp) at 30× coverage completed in ~5 minutes with 32 GB RAM and 8 cores. A moderate genome (~1.3 Gbp) at 30× coverage completed in ~1.5 hours with the same resources. A large genome (~3.3 Gbp) at 30× coverage took several hours.  


PacBio HiFi
- **Script:** `05_variant_calling_hifi.sh`
- **SBATCH wrapper:** `step05_variant_calling_hifi.sbatch`
- **Caller:** `longshot` (v1.0.0)
- **Inputs:** Filtered HiFi BAMs and indexed reference FASTA (`REF.fa` + `REF.fa.fai`).
- **Options:** 
  - `--contigs FILE` to restrict variant calling to a provided list of contigs. Longshot is run per contig and results are concatenated.  
  - `--set-sample` to use BAM basename as the VCF sample ID.  
  - `--outdir` to control output directory (default: `vcf_hifi/`).
- **Outputs:** Per-BAM VCFs in `vcf_hifi/`:  
  `<sample>.hifi.vcf.gz` + index `.tbi` 
- **Effort:** For HiFi, a small genome (~256 Mbp) at 30× coverage completed in ~2.5 hours with 32 GB RAM and 8 cores.



## 6. Build IUPAC diploid consensus fastq

### Generate coverage masks
- **Script:** `06_mask_coverage_outliers.sh`
- **SBATCH wrapper:** `step06_mask_coverage_outliers.sbatch`
- **Purpose:** Generate BED masks of sites with outlier coverage (depth below MIN or above MAX). These masks are later used with `bcftools consensus -m` to mask problematic regions.
- **Inputs:** A plain-text file of BAM paths (filtered BAMs from Step 3b).
- **Outputs:** For each BAM, a BED file in `coverage_masks/` named `<sample>.mask.bed`.
- **Notes:** Default thresholds are `--min 10` and `--max 60`. Adjust as needed (e.g., one-third to twice the mean coverage).  
- **Effort:** Small genomes with 1–5 GB BAMs take 3–5 minutes with ~24 GB RAM on a single CPU. Large genomes with 15 GB BAMs may take up to 7 hours with ~364 GB RAM and 8 threads (though only ~13% CPU used). 3–4 CPUs are generally sufficient.

### Build consensus FASTA
- **Tool:** `bcftools consensus`
- **Purpose:** Apply VCF calls and BED masks to a reference FASTA to build an IUPAC-coded diploid consensus sequence.
- **Inputs:**  
  - Indexed reference FASTA: `-f ref.fa` (+ `ref.fa.fai`)
  - flag to output IUPAC codes `-I`
  - sample name `-s`  
  - BED mask (`-m <sample>.mask.bed`)  
  - Filtered VCF (`.vcf.gz` + `.tbi`)  
- **Outputs:** `-o <sample>.iupac.fa`
- **Effort:** Typically seconds to a minute per sample.

#### Example
`bcftools consensus \
    -f ref.fa \
    -I \
    -s Dryas_octopetala_ERR12370312.hifi.9chrom.singleton \
    -m coverage_masks/Dryas_octopetala_ERR12370312.hifi.9chrom.singleton.mask.bed \
    -o 06_consensus/Dryas_octopetala_ERR12370312.iupac.fa \
    vcf_hifi/Dryas_octopetala_ERR12370312.hifi.vcf.gz` \

### Convert consensus FASTA to FASTQ

- **Tool:** `seqtk`
- **Purpose:** Convert IUPAC consensus FASTA files (`*.iupac.fa`) into FASTQ format with fixed quality scores, which is required as input for PSMC.  
- **Inputs:** `<sample>.iupac.fa` (from Step 6)  
- **Outputs:** `<sample>.iupac.fq`  
- **Notes:** Uses `seqtk seq -F I` where `-F I` sets all bases to quality `I` (Phred Q40).  
- **Effort:** Less than 1 minute per sample.

#### Example
`seqtk seq -F I 06_consensus/Dryas_drummondii_SRR5313982.iupac.fa \
    06_consensus/Dryas_drummondii_SRR5313982.iupac.fq`


## 7. Run PSMC per sample

This stage is done **sample-by-sample**.

#### (Optional) Remove small scaffolds
- **Script:** `remove_small_scaffolds.sh`
- **Purpose:** Length-filter the consensus FASTQ before generating the PSMCFA (keeps reads ≥ MIN bp).
- **Notes:** Requires `seqkit`; writes `<input>.min<MIN>.fq`. :contentReference[oaicite:0]{index=0}

#### Submit PSMC (with bootstraps)
- **Wrapper:** `submit_four_Pfur.sh`
- **Purpose:** Convenience wrapper for Pedicularis furbishiae that submits four parameterizations (Conservative, Standard, Recent-focused, Aggressive-recent) to Slurm.
- **Worker script:** `psmc_boot_multicore.sbatch`
  - **Inputs:** A PSMCFA file and output directory.
  - **Parameters:** Pattern string (e.g., `"4+25*2+4+6"`), mutation rate scaling, iteration settings, number of bootstraps, and number of processors.
  - **Behavior:** Runs the main PSMC analysis and generates bootstrap replicates in parallel, writing results into the specified output directory.
  - **Outputs:**  
    - `<sample>.out_main.psmc` (main run)  
    - `boot/boot_*.psmc` (bootstrap runs)  

#### Plotting PSMC files

`psmc_plot.pl -h`

#### No boots example
`psmc_plot.pl -u 7.0e-9 -R -g 3 -x 100 -Y 2.5 -L -w 4 -T "02_standard_lin" 02_standard_lin 02_standard/out_main.psmc`

#### With boots example
`psmc_plot.pl -u 7.0e-9 -R -g 3 -x 100 -Y 2.5 -L -w 4 -T "02_standard_lin"     02_standard_v2      02_standard/out_main.psmc     02_standard/boot/boot_*.psmc`




---

## Example Workflow for One Sample (Illumina)

```
# 1. Fetch reads from SRA
sbatch --job-name=fetch-SRR5313982 step01_fetch.sbatch SRR5313982

# 2. Clean reads with fastp
sbatch --job-name=fastp-SRR5313982 step02_fastp.sbatch SRR5313982

# 3. Map cleaned reads to reference (Illumina, bwa-mem2)
sbatch --job-name=map-ILL-SRR5313982 \
  map/step03_map_illumina.sbatch \
  Dryas_drummondii_SRR5313982 \
  /path/to/ref.fa \
  /path/to/SRR5313982_R1.clean.fastq.gz \
  /path/to/SRR5313982_R2.clean.fastq.gz \
  illumina
  
## 3b. Filter BAM to keep only primary, mapped, non-dup, non-QC-fail reads
sbatch step03b_filter_bam.sbatch \
  --exclude-flags 3844 \
  map/Dryas_drummondii_SRR5313982.ill.sorted.markdup.bam

# 4. Coverage reporting
bash 04_report_bam_coverage_depth.sh \
  -f bam_list.txt \
  --contigs contigs_of_interest.txt

# 5. Variant calling (Illumina with bcftools)
sbatch step05_variant_calling_ilmn.sbatch \
  -r /path/to/Dryas_drummondii_ref.fa \
  map/Dryas_drummondii_SRR5313982.ill.sorted.markdup.singleton.bam

# 6. Generate coverage mask and build consensus FASTQ
bash 06_mask_coverage_outliers.sh \
  -f bam_list.txt --min 10 --max 60

bcftools consensus \
  -f /path/to/Dryas_drummondii_ref.fa \
  -I \
  -s Dryas_drummondii_SRR5313982.ill.sorted.markdup.singleton \
  -m coverage_masks/Dryas_drummondii_SRR5313982.ill.sorted.markdup.singleton.mask.bed \
  -o 06_consensus/Dryas_drummondii_SRR5313982.iupac.fa \
  vcf_ilmn/Dryas_drummondii_SRR5313982.ilmn.vcf.gz
  
seqtk seq -F I 06_consensus/Dryas_drummondii_SRR5313982.iupac.fa \
  > 06_consensus/Dryas_drummondii_SRR5313982.diploid.fq
  
# 7. Run PSMC per sample

## (Optional) Remove small scaffolds from the diploid FASTQ before PSMCFA generation
bash remove_small_scaffolds.sh \
  06_consensus/Dryas_drummondii_SRR5313982.diploid.fq 500000

## Generate PSMC input
fq2psmcfa -q20 -g10000 \
  06_consensus/Dryas_drummondii_SRR5313982.diploid.min500000.fq \
  > 07_psmc/Dryas_drummondii_SRR5313982.psmcfa

## Submit PSMC runs with bootstraps (example using wrapper)
bash submit_four_Pfur.sh

## Or directly call the worker script for one run
sbatch --cpus-per-task=10 psmc_boot_multicore.sbatch \
  07_psmc/Dryas_drummondii_SRR5313982.psmcfa \
  07_psmc/02_standard "4+25*2+4+6" 10 5 25 100 10

## Plot results without bootstraps
psmc_plot.pl -u 7.0e-9 -R -g 3 -x 100 -Y 2.5 -L -w 4 \
  -T "Dryas drummondii (standard)" \
  07_psmc/02_standard \
  07_psmc/02_standard/out_main.psmc

## Plot results with bootstraps
psmc_plot.pl -u 7.0e-9 -R -g 3 -x 100 -Y 2.5 -L -w 4 \
  -T "Dryas drummondii (standard + boots)" \
  07_psmc/02_standard \
  07_psmc/02_standard/out_main.psmc \
  07_psmc/02_standard/boot/boot_*.psmc
  
```

