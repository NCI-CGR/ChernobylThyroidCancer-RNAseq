# Chernobyl Thyroid Cancer - RNA-seq QC and Alignment
## I. Description
This workflow was used for general QC and alignment of RNA-seq data in the Chernbobyl thyroid cancer study.
Major steps in the workflow are:
1) Trimming of adapters and low-quality reads using trimmomatic
2) Generating QC reports using FASTQC and aggregating results using multiQC
3) Aligning trimmed reads to GRCh38 human reference genome (illumine iGenomes NCBI GRCh38) using STAR
4) Merging reads-count tables of all samples
## II. Dependencies
1) Python
2) Snakemake
3) Trimmomatic
4) Fastqc
5) Multiqc
6) Star
7) R
## III. Input requirements
1) config.yaml
2) sample_names.txt
3) merged fastq files stored in directory: merged_fastq/
4) reference genome sequence and annotation file
