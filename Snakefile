# Assignment 2: Snakefile for Bioinformatics Pipeline
 
rule all:
    input: 
        "results/raw/reference.fasta",
        "results/raw/SRR1972739.fastq",
        "results/qc/SRR1972739_fastqc.zip",
        "results/aligned/aligned.sam",

rule download_reference:
    output:
         reference_fasta = "results/raw/reference.fasta"
    shell:
        """
        echo Downloading reference genome...
        efetch -db nucleotide -id AF086833.2 -format fasta > {output}
        """

rule download_reads:
    output:
        "results/raw/SRR1972739.fastq"
    shell:
        """
        prefetch SRR1972739 -O results/raw/
        fastq-dump -X 10000 --split-3 --outdir results/raw/ results/raw/SRR1972739
        mv results/raw/SRR1972739_1.fastq {output}
        """

rule run_fastqc:
    input:
        "results/raw/SRR1972739.fastq"
    output:
        "results/qc/SRR1972739_fastqc.zip"
    shell:
        """
        echo Running FastQC...
        fastqc -o results/qc {input}
        echo FastQC completed!
        """

rule align_reads:
    input:
        ref="results/raw/reference.fasta",
        reads="results/raw/SRR1972739.fastq"
    output:
        "results/aligned/aligned.sam"
    shell:
        """
        bwa index {input.ref}
        bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:sample1' {input.ref} {input.reads} > {output}
        """


