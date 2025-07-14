rule all:
    input: 
        "results/raw/reference.fasta",
        "results/raw/SRR1972739.fastq",
        "results/qc/SRR1972739_fastqc.zip",
        "results/aligned/aligned.sam",
        "results/aligned/aligned.sorted.bam",
        "results/aligned/validation.txt",
        "results/aligned/dedup.bam",
        "results/aligned/dup_metrics.txt",
        "results/aligned/dedup.bam.bai",
        "results/raw/reference.fasta.fai",
        "results/raw/reference.dict",
        "results/variants/raw_variants.vcf",
        "results/variants/filtered_variants.vcf",
        "results/snpEff/data/reference_db/genes.gbk",
        "results/snpEff/snpEff.config",
       
rule download_reference:
    output:
        "results/raw/reference.fasta"
    shell:
        "efetch -db nucleotide -id AF086833.2 -format fasta > {output}"

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
        "fastqc -o results/qc {input}"

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
rule sam_to_sorted_bam:
    input:
        "results/aligned/aligned.sam"
    output:
        "results/aligned/aligned.sorted.bam"
    shell:
        """
        samtools view -b {input} | samtools sort -o {output}
        """
rule validate_bam:
    input:
        "results/aligned/aligned.sorted.bam"
    output:
        "results/aligned/validation.txt"
    shell:
        """
        gatk ValidateSamFile -I {input} -MODE SUMMARY > {output}
        """
rule mark_duplicates:
    input:
        "results/aligned/aligned.sorted.bam"
    output:
        bam="results/aligned/dedup.bam",
        metrics="results/aligned/dup_metrics.txt"
    shell:
        """
        gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}
        """
rule index_bam:
    input:
        "results/aligned/dedup.bam"
    output:
        "results/aligned/dedup.bam.bai"
    shell:
        "samtools index {input}"     
rule index_reference:
    input:
        "results/raw/reference.fasta"
    output:
        "results/raw/reference.fasta.fai"
    shell:
        "samtools faidx {input}"    

rule create_reference_dict:
    input:
        "results/raw/reference.fasta"
    output:
        "results/raw/reference.dict"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output}"

rule call_variants:
    input:
        bam="results/aligned/dedup.bam",
        bai="results/aligned/dedup.bam.bai",
        ref="results/raw/reference.fasta"
    output:
        "results/variants/raw_variants.vcf"
    shell:
        """
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output}
        """        
rule filter_variants:
    input:
        vcf="results/variants/raw_variants.vcf",
        ref="results/raw/reference.fasta"
    output:
        "results/variants/filtered_variants.vcf"
    shell:
        """
        gatk VariantFiltration \
          -R {input.ref} \
          -V {input.vcf} \
          -O {output} \
          --filter-expression "QD < 2.0 || FS > 60.0" \
          --filter-name FILTER
        """

rule download_genbank:
    input:
        "results/raw/reference.fasta"
    output:
        "results/snpEff/data/reference_db/genes.gbk"
    shell:
        "efetch -db nucleotide -id AF086833.2 -format genbank > {output}"

rule build_snpeff_config:
    input:
        fasta="results/raw/reference.fasta",
        genbank="results/snpEff/data/reference_db/genes.gbk"
    output:
        "results/snpEff/snpEff.config"
    run:
        fasta = os.path.abspath(input.fasta)
        genbank = os.path.abspath(input.genbank)
        with open(output[0], "w") as f:
            f.write(f"""reference_db.genome : reference_db
reference_db.fa : {fasta}
reference_db.genbank : {genbank}
""")
