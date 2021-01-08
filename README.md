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
* [merged fastq files stored in directory: merged_fastq/](https://github.com/NCI-CGR/ChernobylThyroidCancer-RNAseq/tree/main/merged_fastq)
* reference genome sequence and annotation files
## IV. Working directory structure
```bash
.
├── config.yaml
├── merged_fastq
│   ├── {sample}_merged_R1.fastq.gz
│   └── {sample}_merged_R2.fastq.gz
├── merge.R
├── posttrim_qc
│   ├── posttrim_qc_multiqc_report.html
│   ├── {sample}_filtered_1P_fastqc.html
│   ├── {sample}_filtered_1P_fastqc.zip
│   ├── {sample}_filtered_2P_fastqc.html
│   └── {sample}_filtered_2P_fastqc.zip
├── pretrim_qc
│   ├── pretrim_qc_multiqc_report.html
│   ├── {sample}_merged_R1_fastqc.html
│   ├── {sample}_merged_R1_fastqc.zip
│   ├── {sample}_merged_R2_fastqc.html
│   └── {sample}_merged_R2_fastqc.zip
├── reads_count
│   └── reads_count.csv
├── run.sh
├── sample_names.txt
├── Snakefile
├── snakemake.batch
├── star_align
│   ├── log
│   │   └── star_align_multiqc_report.html
│   └── {sample}
│       ├── {sample}Aligned.sortedByCoord.out.bam
│       ├── {sample}Log.final.out
│       └── {sample}ReadsPerGene.out.tab
├── star_index
│   ├── chrLength.txt
│   ├── chrNameLength.txt 
│   ├── chrName.txt 
│   ├── chrStart.txt 
│   ├── exonGeTrInfo.tab 
│   ├── exonInfo.tab 
│   ├── geneInfo.tab 
│   ├── Genome 
│   ├── genomeParameters.txt 
│   ├── SA 
│   ├── SAindex 
│   ├── sjdbInfo.txt 
│   ├── sjdbList.fromGTF.out.tab 
│   ├── sjdbList.out.tab 
│   └── transcriptInfo.tab 
└── trimmed
    ├── {sample}_filtered_1P.fq.gz
    ├── {sample}_filtered_1U.fq.gz
    ├── {sample}_filtered_2P.fq.gz
    └── {sample}_filtered_2U.fq.gz
```
