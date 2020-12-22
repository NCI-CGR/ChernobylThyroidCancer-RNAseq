# Chernobyl Thyroid Cancer - RNA-seq QC and Alignment
## I. Description
This workflow was used for general QC and alignment of RNA-seq data in the Chernbobyl thyroid cancer study.
Major steps in the workflow are:
1) Trimming of adapters and low-quality reads using trimmomatic/0.36
2) Generating QC reports using FASTQC/0.11.5 and aggregating results using multiQC
3) Aligning trimmed reads to GRCh38 human reference genome (illumine iGenomes NCBI GRCh38) using STAR/2.5.4a
4) Merging reads-count tables of all samples
