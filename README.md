# Chernobyl Thyroid Cancer - RNA-seq QC and Alignment
## I. Description
This workflow was used for general QC and alignment of RNA-seq data in the Chernbobyl thyroid cancer study.
Major steps in the workflow are:
1) Trimming of adapters and low-quality reads using trimmomatic
2) Generating QC reports using FASTQC and aggregating results using multiQC
3) Aligning trimmed reads to GRCh38 human reference genome (illumine iGenomes NCBI GRCh38) using STAR
4) Merging reads-count tables of all samples
## II. Dependencies
* [Python](https://www.python.org)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Multiqc](https://multiqc.info)
* [Star](https://github.com/alexdobin/STAR)
* [R](https://www.r-project.org)
## III. Input requirements
* [config.yaml](https://github.com/NCI-CGR/ChernobylThyroidCancer-RNAseq/blob/main/config.yaml)
* [sample_names.txt](https://github.com/NCI-CGR/ChernobylThyroidCancer-RNAseq/blob/main/sample_names.txt)
* merged fastq files stored in directory: merged_fastq/
* reference genome sequence and annotation files
