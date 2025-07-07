# Assignment 2: Snakefile for Bioinformatics Pipeline
 
rule all:
    input: 
        "results/raw/reference.fasta",
    
 
rule download_reference:
    output:
         reference_fasta = "results/raw/reference.fasta"
    shell:
        """
        echo Downloading reference genome...
        efetch -db nucleotide -id AF086833.2 -format fasta > {output}
        """