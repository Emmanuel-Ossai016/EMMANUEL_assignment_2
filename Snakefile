# Assignment 2: Snakefile for Bioinformatics Pipeline
 
rule all:
    input: 
        "results/raw/reference.fasta",
        "results/raw/SRR1972739.fastq",
    
 
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