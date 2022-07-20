## snakemake pipeline for ChIP-seq analysis
# Genome files 
genome_index="/home/user/genomes/bowtie2Index/mm10/DNA/genome"
genome_gtf="/home/user/genomes/mm10.refGene.gtf"
genome_gtf_gene="/home/user/genomes/mm10.refGene.transcript.gtf"
# input files
species="Mouse/Mus musculus"
reads=["_R1","_R2"]    #paired-end sequencing
ext=".fastq.gz"
path_origin="/media/temp/origindata"    #path to raw sequencing data

sample_input=["R21072088-P6-BMI1-input-P6-BMI1-input_combined",
            "R21088164-H3K27me3input-H3K27me3input_combined"]

sample_chip=["R21007067-P5-B-P5-B_combined",
            "R19054399-libmix20191202-2-bmil-ChIP_combined",
            "R21042415-p6-H2Aub-p6-H2Aub_combined",
            "H2Aub_S20_L003",
            "R21042416-p6-H3k27me3-p6H3k27me3-1_combined",
            "H3K27me3_S21_L003",
            "R21094462-K27ac-K27ac_combined",
            "K27ac-c3_S0_L000",
            "R21094463-KD-H2Aub-KD-H2Aub_combined",
            "KD-ab_S15_L004",]
##
import os
import sys
from pathlib import Path

samples=sample_chip+sample_input
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
    pass
print("#"*10)
print("Species is:",species)
print("All samples:", samples)
print("#"*10)
##
rule all:
    input:
        expand("fastQC/{sample}{read}_fastqc.html",sample=samples,read=reads),
        expand("alignment/{sample}.sorted.bam",sample=samples),
        expand("alignment/{sample}.sorted.bam.bai",sample=samples),
        expand("deeptools/bigwig/{sample}.filter.dedup.bw", sample=samples),
        expand("deeptools/fragsize/{sample}_fragsize.png", sample=samples),
        expand("deeptools/profile/{sample}_gene.bed",sample=samples),
        expand("deeptools/profile/{sample}_gene.png",sample=samples),
        "deeptools/plotfingerprint.png",
        "deeptools/multibamsummary.counts.txt",
        "deeptools/correlation.png"

rule trimming:
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

rule fastQC:
    input:
        "trimdata/{sample}{read}"+ext
    output:
        "fastQC/{sample}{read}_fastqc.html"
    threads: 24
    log:
        "fastQC/logs/{sample}{read}.log"
    shell:
        "fastqc -o fastQC -t {threads} -q {input} 2> {log}"

rule bowtie2:
    input:
        r1="trimdata/{sample}"+reads[0]+ext,
        r2="trimdata/{sample}"+reads[1]+ext
    output:
        bam="alignment/{sample}.sorted.bam",
        summary="alignment/{sample}.summary"
    threads: 24
    params:
        index=genome_index
    log:
        "alignment/logs/{sample}.align.log"
    shell:
        "bowtie2 -q -X 1000 -x {params.index} -1 {input.r1} -2 {input.r2} -p {threads} 2> {output.summary} "
        "| samtools view -Sbh - | samtools sort -@ {threads} -o {output.bam} 2> {log}"

rule samindex:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        "alignment/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule filter:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        "alignment/filter/{sample}.filter.bam"
    threads: 24
    log:
        "alignment/filter/logs/{sample}.filter.log"
    shell:
        "samtools view -b -F 1804 -f 2 -q 30 -@ {threads} {input} > {output} 2> {log}"

rule deduplicates:
    input:
        "alignment/filter/{sample}.filter.bam"
    output:
        bam="alignment/filter/{sample}.filter.dedup.bam",
        matrix="alignment/filter/{sample}.filter.dedup.txt"
    threads: 16
    log:
        "alignment/filter/logs/{sample}.dedup.log"
    shell:
        "picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.matrix} "
        "REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT 2> {log}"
    
rule fragmentszie:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        png="deeptools/fragsize/{sample}_fragsize.png",
        hist="deeptools/fragsize/{sample}_fragsize.hist.txt"
    threads: 24
    log:
        "deeptools/logs/{sample}.fragsize.log"
    shell:
        "bamPEFragmentSize -b {input} -hist {output.png} --outRawFragmentLengths {output.hist} --maxFragmentLength 1000"

rule bamcoverage:
    input:
        "alignment/filter/{sample}.filter.dedup.bam"
    output:
        "deeptools/bigwig/{sample}.filter.dedup.bw"
    params:
        black=blacklist
    threads: 16
    log:
        "deeptools/logs/{sample}.bamcoverage.log"
    shell:
        "bamCoverage -b {input} -o {output} -p {threads} "
        "--binSize 10 --smoothLength 30 --normalizeUsing RPKM --ignoreForNormalization chrX chrY chrM --extendReads 2> {log}"

rule profile:
    input:
        "deeptools/bigwig/{sample}.filter.dedup.bw"
    output:
        bed="deeptools/profile/{sample}_gene.bed",
        matrix="deeptools/profile/{sample}_gene.matrix.gz"
    params:
        ref=genome_gtf_gene
    threads: 24
    log:
        "deeptools/logs/{sample}.profile.log"
    shell:
        "computeMatrix reference-point --referencePoint TSS -b 2500 -a 2500 -R {params.ref} "
        "-S {input} -o {output.matrix} --outFileSortedRegions {output.bed} --skipZeros -p {threads} 2> {log}"

rule plotprofile:
    input:
        "deeptools/profile/{sample}_gene.matrix.gz"
    output:
        "deeptools/profile/{sample}_gene.png"
    shell:
        "plotHeatmap -m {input} -out {output} --dpi 360"

rule plotfingerprint:
    input:
        expand("alignment/filter/{sample}.filter.dedup.bam",sample=samples)
    output:
        matrix="deeptools/plotfingerprint.txt",
        png="deeptools/plotfingerprint.png"
    threads: 24
    log:
        "deeptools/logs/plotfingerprint.log"
    shell:
        "plotFingerprint -b {input} -o {output.png} --outQualityMetrics {output.matrix} "
        "-p {threads} --ignoreDuplicates --extendReads 2> {log}"

rule correlation:
    input:
        expand("alignment/filter/{sample}.filter.dedup.bam",sample=samples)
    output:
        counts="deeptools/multibamsummary.counts.txt",
        npz="deeptools/multibamsummary.counts.npz"
    threads: 24
    log:
        "deeptools/logs/multibamsummary.log"
    shell:
        "multiBamSummary bins --bamfiles {input} --outRawCounts {output.counts} -o {output.npz} --extendReads --binSize 2000 -p {threads} 2> {log}; "

rule plotcorrelation:
    input:
        "deeptools/multibamsummary.counts.npz"
    output:
        "deeptools/correlation.png"
    shell:
        "plotCorrelation -in {input} -o {output} --corMethod pearson --whatToPlot scatterplot --skipZeros --removeOutliers --log1p"
#...
