SAMPLES = ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6']

rule star_align:
    input:
        R1="data/{sample}_R1.fastq.gz",
        R2="data/{sample}_R2.fastq.gz",
        index="ref/"
    output:
        bam="alignment/{sample}.Transcriptome.out.bam"
    log:
        "logs/star/{sample}.log"
    conda:
        "envs/star.yaml"
    shell:
        """
        STAR --runThreadN 8 --genomeDir {input.index} --readFilesIn {input.R1} {input.R2} --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat > {log} 2>&1
        mv alignment/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}
        """

rule salmon_quant:
    input:
        bam="alignment/{sample}.Transcriptome.out.bam"
    output:
        dir="quant/{sample}_quant/"
    log:
        "logs/salmon/{sample}.log"
    conda:
        "envs/salmon.yaml"
    shell:
        """
        salmon quant -t homo_index.fa -l A -a {input.bam} -o {output.dir} > {log} 2>&1
        """

rule all:
    input:
        expand("alignment/{sample}.Transcriptome.out.bam", sample=SAMPLES),
        expand("quant/{sample}_quant", sample=SAMPLES)
