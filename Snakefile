### This pipeline is for general QC and alignment of RNA-seq data

## vim: ft=python
import sys
import os
import glob
import itertools

shell.prefix("set -eo pipefail; ")
localrules: all

# define wildcards
def parse_sampleID(fname):
    return fname.split('/')[-1].split('_')[0]

file = sorted(glob.glob('merged_fastq/*.fastq.gz'), key=parse_sampleID)

d = {}
for key, value in itertools.groupby(file, parse_sampleID):
    d[key] = list(value)

rule all:
    input:
          expand("star_align/{sample}/{sample}Aligned.sortedByCoord.out.bam",sample=d.keys()),
          "pretrim_qc/preQC_multiqc_report.html",
          "posttrim_qc/postQC_multiqc_report.html",
          "star_align/log/star_align_multiqc_report.html",
          "reads_count/reads_count.csv"

rule trimmomatic:
    input:
          "merged_fastq/{sample}_merged_R1.fastq.gz",
          "merged_fastq/{sample}_merged_R2.fastq.gz"
    output:
          "trimmed/{sample}_filtered_1P.fq.gz",
          "trimmed/{sample}_filtered_1U.fq.gz",
          "trimmed/{sample}_filtered_2P.fq.gz",
          "trimmed/{sample}_filtered_2U.fq.gz"
    threads: 16 
    shell:
          """
          java -classpath /usr/local/apps/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 -threads 24 {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>log/{wildcards.sample}_trim.err 
          """

rule pretrim_qc:
    input:
          "merged_fastq/{sample}_merged_R1.fastq.gz",
          "merged_fastq/{sample}_merged_R2.fastq.gz"
    output: 
          "pretrim_qc/{sample}_merged_R1_fastqc.zip",
          "pretrim_qc/{sample}_merged_R1_fastqc.html",
          "pretrim_qc/{sample}_merged_R2_fastqc.zip",
          "pretrim_qc/{sample}_merged_R2_fastqc.html"
    threads: 10
    shell:
          """
          fastqc {input} -o pretrim_qc -f fastq --noextract 2>log/{wildcards.sample}_preqc.err
          """

rule posttrim_qc:
    input:
          "trimmed/{sample}_filtered_1P.fq.gz",
          "trimmed/{sample}_filtered_2P.fq.gz"
    output:
          "posttrim_qc/{sample}_filtered_1P_fastqc.zip",
          "posttrim_qc/{sample}_filtered_1P_fastqc.html",
          "posttrim_qc/{sample}_filtered_2P_fastqc.zip",
          "posttrim_qc/{sample}_filtered_2P_fastqc.html"
    threads: 10
    shell:
          """
          fastqc {input} -o posttrim_qc -f fastq --noextract 2>log/{wildcards.sample}_postqc.err
          """

rule star_index:
    input:
          "/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa",
          "/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
    output:
          "star_index/complete.txt"
    threads: 24
    shell:
          """
          STAR --runThreadN 24 --runMode genomeGenerate --genomeDir star_index --sjdbGTFfile {input[1]} --sjdbOverhang 149 --genomeFastaFiles {input[0]} 2>log/star_index.err
          touch {output} 
          """

rule star_align:
    input:
          "/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf",
          "trimmed/{sample}_filtered_1P.fq.gz",
          "trimmed/{sample}_filtered_2P.fq.gz",
          "star_index/complete.txt"
    output:
          "star_align/{sample}/{sample}Aligned.sortedByCoord.out.bam",
          "star_align/{sample}/{sample}ReadsPerGene.out.tab",
          "star_align/{sample}/{sample}Log.final.out"
    threads: 24
    params:
          index="star_index"
    shell:
          """
          STAR --runThreadN 24 --genomeDir {params.index} --readFilesIn {input[1]} {input[2]} --outFileNamePrefix star_align/{wildcards.sample}/{wildcards.sample} --readFilesCommand zcat --sjdbGTFfile {input[0]} --sjdbOverhang 149 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --twopassMode Basic 2>log/{wildcards.sample}_star_align.err
          """

rule multiqc:
    input:
          expand("pretrim_qc/{sample}_merged_R1_fastqc.html",sample=d.keys()),
          expand("posttrim_qc/{sample}_filtered_1P_fastqc.html",sample=d.keys()),
          expand("star_align/{sample}/{sample}Log.final.out",sample=d.keys())
    output:
          "pretrim_qc/preQC_multiqc_report.html",
          "posttrim_qc/postQC_multiqc_report.html",
          "star_align/log/star_align_multiqc_report.html"
    threads: 8
    shell:
          """
          multiqc pretrim_qc/. --title preQC -o pretrim_qc 2>log/multiqc_preqc.err
          multiqc posttrim_qc/. --title postQC -o posttrim_qc/ 2>log/multiqc_postqc.err
          mkdir star_align/log
          cp star_align/*/*Log.final.out star_align/log
          multiqc star_align/log/. --title star_align -o star_align/log 2>log/multiqc_star.err
          """

##library is reverse stranded
rule merge:
    input:  
          expand("star_align/{sample}/{sample}ReadsPerGene.out.tab",sample=d.keys())
    output:
          "reads_count/reads_count.csv"
    shell:
          """
          Rscript merge.R 2>log/merge_count.err
          """

