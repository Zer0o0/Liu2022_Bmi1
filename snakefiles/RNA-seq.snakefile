##snakemake pipeline for RNA-seq analysis
# Genome files
genome_index="/home/user/genomes/hisat2Index/mm10/genome"
genome_gtf="/home/user/genomes/mm10.refGene.gtf"
# input files
species="Mouse/Mus musculus"
reads=["_R1","_R2"]
ext=".fastq.gz"
path_origin="/media/temp/origindata"

samples=["KD_S25_L001",
        "WT_S0_L000",
        "R21070857-P6-WT-C1-P6-WT-C1_combined",
        "R21070858-P6-KD-48H-C3-P6-KD-48H-C3_combined",
        "R21042411-P6-WT-P6-WT-C2_combined",
        "R21042412-P6-KD-48H-P6-KD-48H_combined",]
##
import os
import sys
from pathlib import Path

p = Path(path_origin)
ss=list(p.glob("**/*"+ext))
l=len(ext)+len(reads[0])
try:
    Path("rawdata").mkdir(parents=True, exist_ok=True)
    pwd=os.getcwd()
    for s in ss:
        f=s.name
        if f[0:-l] in samples:
            os.symlink(s,os.path.join(pwd,"rawdata",f))
except:
    print("Something wrong in pre-pipeline!!!")
    #sys.exit(1)
    pass
print("#"*10)
print("Species is:")
print("All samples:", samples)
print("#"*10)
##
rule all:
    input:
        expand("fastQC/{sample}{read}_fastqc.html",sample=samples,read=reads),
        expand("alignment/{sample}.filter.sorted.bam.bai",sample=samples),
        expand("featurecounts/{sample}.count.txt",sample=samples),
        expand("cufflinks/{sample}/genes.fpkm_tracking",sample=samples)

rule trimgalore:
    input:
        r1="rawdata/{sample}"+reads[0]+ext,
        r2="rawdata/{sample}"+reads[1]+ext
    output:
        r1="trimdata/{sample}"+reads[0]+ext,
        r2="trimdata/{sample}"+reads[1]+ext
    params:
        tmp1="trimdata/{sample}"+reads[0]+"_val_1.fq.gz",
        tmp2="trimdata/{sample}"+reads[1]+"_val_2.fq.gz"
    threads: 4
    log:
        "trimdata/logs/{sample}.log"
    shell:
        "trim_galore --paired --stringency 3 -q 20 --length 20 -o trimdata --cores {threads} {input.r1} {input.r2} 2> {log}; "
        "mv {params.tmp1} {output.r1}; "
        "mv {params.tmp2} {output.r2}"

rule fastqc:
    input:
        "trimdata/{sample}{read}"+ext
    output:
        "fastQC/{sample}{read}_fastqc.html"
    log:
        "fastQC/logs/{sample}{read}.log"
    shell:
        "fastqc -o fastQC -t {threads} -q {input} 2> {log}"

rule hisat2:
    input:
        r1="trimdata/{sample}"+reads[0]+ext,
        r2="trimdata/{sample}"+reads[1]+ext
    output:
        bam="alignment/{sample}.bam",
        summary="alignment/{sample}.summary"
    threads: 24
    params:
        index=genome_index,
        splice="alignment/{sample}.splice.txt",
        met="alignment/{sample}.matrix.txt"
    log:
        "alignment/logs/{sample}.align.log"
    shell:
        "hisat2 -q --rna-strandness RF --dta-cufflinks -x {params.index} -1 {input.r1} -2 {input.r2} -p {threads} --novel-splicesite-outfile {params.splice} --summary-file {output.summary} --met-file {params.met} | "
        "samtools view -Sbh - > {output.bam} 2> {log}"

rule filter:
    input:
        "alignment/{sample}.bam"
    output:
        "alignment/{sample}.filter.sorted.bam"
    threads: 16
    log:
        "alignment/logs/{sample}.filter.log"
    shell:
        "samtools view -b -F 780 -f 2 -q 30 -@ {threads} {input} | "
        "samtools sort -@ {threads} -o {output} 2> {log}"

rule bamindex:
    input:
        "alignment/{sample}.filter.sorted.bam"
    output:
        "alignment/{sample}.filter.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule featurecounts:
    input:
        "alignment/{sample}.filter.sorted.bam"
    output:
        "featurecounts/{sample}.count.txt"
    params:
        gtf=genome_gtf
    log:
        "featurecounts/logs/{sample}.count.log"
    shell:
        "featureCounts -s 2 -p -B -t exon -g gene_id -a {params.gtf} -o {output} {input} 2> {log}"

rule cufflink:
    input:
        "alignment/{sample}.filter.sorted.bam"
    output:
        "cufflinks/{sample}/genes.fpkm_tracking"
    threads: 16
    params:
        gtf=genome_gtf,
        out="cufflinks/{sample}"
    log:
        "cufflinks/logs/{sample}.log"
    shell:
        "cufflinks -p {threads} -o {params.out} -g {params.gtf} --library-type fr-firststrand {input} 2> {log}"
