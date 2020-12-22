## vim: ft=python
import sys
import os
from snakemake.utils import R

# 1) adapter trim
# 2) QC
# 3) star index
# 4) star alignment and quantification

shell.prefix("set -eo pipefail; ")
configfile:"config.yaml"
localrules: all

report:"report/workflow.rst"

rule all:
   input:
         expand("/data/DCEG_Chernobyl/phase3/totalRNA/trimmed/{sample}_filtered_1P.fq.gz",sample=config["samples"]),
         expand("/data/DCEG_Chernobyl/phase3/totalRNA/trimmed/{sample}_filtered_2P.fq.gz",sample=config["samples"]),
         expand("/data/DCEG_Chernobyl/phase3/totalRNA/pretrim_qc/{sample}_merged_R1_fastqc.zip",sample=config["samples"]),
         expand("/data/DCEG_Chernobyl/phase3/totalRNA/pretrim_qc/{sample}_merged_R2_fastqc.zip",sample=config["samples"]),
         expand("/data/DCEG_Chernobyl/phase3/totalRNA/posttrim_qc/{sample}_filtered_1P_fastqc.zip",sample=config["samples"]),
         expand("/data/DCEG_Chernobyl/phase3/totalRNA/posttrim_qc/{sample}_filtered_2P_fastqc.zip",sample=config["samples"]),
         "/data/daij/project/RD168_Chernobyl_TN-Pairs_miRNA_HiSeq/total_RNA/star_index/geneInfo.tab",
         expand("/data/DCEG_Chernobyl/phase3/totalRNA/star_align/{sample}/{sample}Aligned.sortedByCoord.out.bam",sample=config["samples"]),
         expand("/data/DCEG_Chernobyl/phase3/totalRNA/star_align/{sample}/{sample}ReadsPerGene.out.tab",sample=config["samples"]),
         "/data/DCEG_Chernobyl/phase3/totalRNA/pretrim_qc/phase3_totalRNA_preQC_multiqc_report.html",
          "/data/DCEG_Chernobyl/phase3/totalRNA/posttrim_qc/phase3_totalRNA_postQC_multiqc_report.html",
          "/data/DCEG_Chernobyl/phase3/totalRNA/star_align/log/phase3_totalRNA_star_align_multiqc_report.html",
          "/data/DCEG_Chernobyl/phase3/totalRNA/reads_count/phase3_totalRNA_reads_count.csv"


rule trimmomatic:
   input:
     "/data/DCEG_Chernobyl/phase3/totalRNA/merged_fastq/{sample}_merged_R1.fastq.gz",
     "/data/DCEG_Chernobyl/phase3/totalRNA/merged_fastq/{sample}_merged_R2.fastq.gz"
   output:
     "/data/DCEG_Chernobyl/phase3/totalRNA/trimmed/{sample}_filtered_1P.fq.gz",
     "/data/DCEG_Chernobyl/phase3/totalRNA/trimmed/{sample}_filtered_1U.fq.gz",
     "/data/DCEG_Chernobyl/phase3/totalRNA/trimmed/{sample}_filtered_2P.fq.gz",
     "/data/DCEG_Chernobyl/phase3/totalRNA/trimmed/{sample}_filtered_2U.fq.gz"

   threads:16 
   shell:
         """
         module load trimmomatic
         java -classpath /usr/local/apps/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 -threads 24 {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
         """

rule pretrim_qc:
    input:
     "/data/DCEG_Chernobyl/phase3/totalRNA/merged_fastq/{sample}_merged_R1.fastq.gz",
     "/data/DCEG_Chernobyl/phase3/totalRNA/merged_fastq/{sample}_merged_R2.fastq.gz"
    output: 
          "/data/DCEG_Chernobyl/phase3/totalRNA/pretrim_qc/{sample}_merged_R1_fastqc.zip",
          "/data/DCEG_Chernobyl/phase3/totalRNA/pretrim_qc/{sample}_merged_R1_fastqc.html",
          "/data/DCEG_Chernobyl/phase3/totalRNA/pretrim_qc/{sample}_merged_R2_fastqc.zip",
          "/data/DCEG_Chernobyl/phase3/totalRNA/pretrim_qc/{sample}_merged_R2_fastqc.html"

    threads: 10
    shell:
         """
         module load fastqc
         fastqc {input} -o /data/DCEG_Chernobyl/phase3/totalRNA/pretrim_qc -f fastq --noextract
         """

rule posttrim_qc:
    input:
     "/data/DCEG_Chernobyl/phase3/totalRNA/trimmed/{sample}_filtered_1P.fq.gz",
     "/data/DCEG_Chernobyl/phase3/totalRNA/trimmed/{sample}_filtered_2P.fq.gz"
    output:
          "/data/DCEG_Chernobyl/phase3/totalRNA/posttrim_qc/{sample}_filtered_1P_fastqc.zip",
          "/data/DCEG_Chernobyl/phase3/totalRNA/posttrim_qc/{sample}_filtered_1P_fastqc.html",
          "/data/DCEG_Chernobyl/phase3/totalRNA/posttrim_qc/{sample}_filtered_2P_fastqc.zip",
          "/data/DCEG_Chernobyl/phase3/totalRNA/posttrim_qc/{sample}_filtered_2P_fastqc.html"

    threads: 10
    shell:
         """
         module load fastqc
         fastqc {input} -o /data/DCEG_Chernobyl/phase3/totalRNA/posttrim_qc -f fastq --noextract
         """

rule star_index:
    input:
        "/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa",
        "/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
    output:
        "/data/daij/project/RD168_Chernobyl_TN-Pairs_miRNA_HiSeq/total_RNA/star_index/geneInfo.tab"
    threads: 24
    shell:
        """
        module load STAR
        STAR --runThreadN 24 --runMode genomeGenerate --genomeDir /data/daij/project/RD168_Chernobyl_TN-Pairs_miRNA_HiSeq/total_RNA/star_index --sjdbGTFfile {input[1]} --sjdbOverhang 149 --genomeFastaFiles {input[0]} 
        """

rule star_align:
        input:
           "/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf",
           "/data/DCEG_Chernobyl/phase3/totalRNA/trimmed/{sample}_filtered_1P.fq.gz",
           "/data/DCEG_Chernobyl/phase3/totalRNA/trimmed/{sample}_filtered_2P.fq.gz",
        output:
           "/data/DCEG_Chernobyl/phase3/totalRNA/star_align/{sample}/{sample}Aligned.sortedByCoord.out.bam",
           "/data/DCEG_Chernobyl/phase3/totalRNA/star_align/{sample}/{sample}ReadsPerGene.out.tab"
        threads:24
        params:
           index="/data/daij/project/RD168_Chernobyl_TN-Pairs_miRNA_HiSeq/total_RNA/star_index"
        shell:
           """
           cd /data/DCEG_Chernobyl/phase3/totalRNA/star_align/
           module load STAR
           STAR --runThreadN 24 --genomeDir /data/daij/project/RD168_Chernobyl_TN-Pairs_miRNA_HiSeq/total_RNA/star_index --readFilesIn {input[1]} {input[2]} --outFileNamePrefix /data/DCEG_Chernobyl/phase3/totalRNA/star_align/{wildcards.sample}/{wildcards.sample} --readFilesCommand zcat --sjdbGTFfile {input[0]} --sjdbOverhang 124 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --twopassMode Basic
           """

rule multiqc:
       input:
          expand("/data/DCEG_Chernobyl/phase3/totalRNA/pretrim_qc/{sample}_merged_R1_fastqc.html",sample=config["samples"]),
          expand("/data/DCEG_Chernobyl/phase3/totalRNA/posttrim_qc/{sample}_filtered_1P_fastqc.html",sample=config["samples"]),
          expand("/data/DCEG_Chernobyl/phase3/totalRNA/star_align/{sample}/{sample}Log.final.out",sample=config["samples"])
       output:
          "/data/DCEG_Chernobyl/phase3/totalRNA/pretrim_qc/phase3_totalRNA_preQC_multiqc_report.html",
          "/data/DCEG_Chernobyl/phase3/totalRNA/posttrim_qc/phase3_totalRNA_postQC_multiqc_report.html",
          "/data/DCEG_Chernobyl/phase3/totalRNA/star_align/log/phase3_totalRNA_star_align_multiqc_report.html"
       threads:8
       shell:
          """
          cd /data/DCEG_Chernobyl/phase3/totalRNA/pretrim_qc/
          module load multiqc
          multiqc . --title phase3_totalRNA_preQC
          cd /data/DCEG_Chernobyl/phase3/totalRNA/posttrim_qc/
          multiqc . --title phase3_totalRNA_postQC
          cd /data/DCEG_Chernobyl/phase3/totalRNA/star_align
          mkdir log
          cp {wildcards.sample}/{wildcards.sample}Log.final.out ./log
          cd log
          multiqc . --title phase3_totalRNA_star_align
          """

##library is first stranded
rule merge:
   input:  expand("/data/DCEG_Chernobyl/phase3/totalRNA/star_align/{sample}/{sample}ReadsPerGene.out.tab",sample = config["samples"])
   output:"/data/DCEG_Chernobyl/phase3/totalRNA/reads_count/phase3_totalRNA_reads_count.csv"
   shell:
          """
          Rscript merge.R
          """

