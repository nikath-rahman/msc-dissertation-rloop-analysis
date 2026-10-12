# msc-dissertation-rloop-analysis

# R-loops at DNA double-strand breaks: analysis scripts

Scripts for the MSc Bioinformatics and Genomics dissertation
"Investigating the role of R-loops in the pathogenesis of neurodevelopmental disorders", Cardiff University, 2026.

Two datasets were analysed against the T2T-CHM13v2.0 reference genome:

- **spKAS-seq** (DIvA cells): one enriched library (Enrich_S2) and one input
  library (Input_S1).
- **RIAN-seq**: 13 libraries from DIvA and RPE1 cells, covering UPF1 depletion
  or knockout, DSB induction (tamoxifen, etoposide) and RNase H treatment.

Shell scripts ran as SLURM jobs on the Falcon cluster (Supercomputing Wales);
R scripts ran locally. Each script takes the project root from the `DISS`
environment variable, for example:

```bash
export DISS=/path/to/dissertation
sbatch 01_qc_trimming/01_fastqc_spkas.sh
```

## Order of use

| Folder | Scripts | Purpose |
| --- | --- | --- |
| `01_qc_trimming` | 01-04 | FastQC and Trim Galore for both datasets (`--nextera` for the Tn5-tagmented RIAN-seq libraries) |
| `02_alignment` | 01 | Bowtie2 index build and alignment to CHM13v2.0 |
| `03_filtering` | 01-03 | MAPQ >= 20 and proper-pair filtering; fixmate and duplicate removal for RIAN-seq; insert-size summary |
| `04_peak_calling` | 01-02 | MACS3: spKAS-seq against its input control; RIAN-seq without a control, with per-sample `--extsize` |
| `05_dsb_proximity` | 01-04 | DSB BED file; overlap and 1 kb window tests for spKAS-seq and DIvA; genome-wide summary for RPE1 |
| `06_coverage_tracks` | 01-02 | deepTools bamCoverage bigWig tracks (RPGC for spKAS-seq, CPM for RIAN-seq) |
| `07_R_analysis` | see folder | regioneR permutation tests, ChIPseeker annotation and figures |

## Software

FastQC 0.12.1, Trim Galore 0.6.10, cutadapt 4.9, Bowtie2 2.5.4,
SAMtools 1.21, MACS3 3.0.4, BEDTools 2.31.1, deepTools [version],
R version 4.5.1 with regioneR 1.42.0 and ChIPseeker 1.46.1.
DSB coordinates were converted from hg19 to CHM13v2.0 with UCSC liftOver.

## Data availability

Raw sequencing data (FASTQ) and alignments (BAM) are not included. They were
supplied by the School of Biosciences, Cardiff University. 
