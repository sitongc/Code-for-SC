SAMPLES = ['sample1', 'sample2', 'sample3']

rule star_align:
    input:
        R1="data/{sample}_R1.fastq.gz",
        R2="data/{sample}_R2.fastq.gz",
        index="/home/ec2-user"
    output:
        bam="alignment/{sample}.sortedByCoord.out.bam"
    log:
        "logs/star/{sample}.log"
    conda:
        "envs/star.yaml"
    shell:
        """
    STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --runThreadN 8 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 --alignIntronMax 299999 --genomeDir {input.index}  --sjdbGTFfile gencode.v42lift37.annotation.gtf --outFileNamePrefix ./ --readFilesIn  {input.R1} {input.R2}
        mv alignment/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}
        """


rule prepare_darts_directory:
    output:
        directory="{sample}_darts"
    shell:
        "mkdir -p {output.directory}"

rule rmats_count:
    input:
        bam="{sample}.sortedByCoord.out.bam",
        b1="b1.txt",
        b2="b2.txt",
        gtf="ref/gencode.v42.chr_patch_hapl_scaff.annotation.gtf",
        sample_dir="{sample}_darts"
    output:
        touch("{sample_dir}/rmats_count_complete.txt")
    params:
        nthread=24,
        readLength=101
    shell:
        """
        echo {input.bam} > {input.b1}
        Darts_BHT rmats_count --b1 {input.b1} --b2 {input.b2} --gtf {input.gtf} --od {output.sample_dir} -t paired --nthread {params.nthread} --readLength {params.readLength}
        """

rule bayes_infer_A5SS:
    input:
        rmats_count_complete="{sample}_darts/rmats_count_complete.txt",
        rmats_count_output="JC.raw.input.A5SS.txt"
    output:
        "A5SS_BHT.results.xlsx"
    params:
        annot="fromGTF.A5SS.txt",
        thread=8,
        c=0.05
    shell:
        """
        rm b1.txt
        cd {sample}_darts/
        Darts_BHT bayes_infer --rmats-count {input.rmats_count_output} --od ./ --annot {params.annot} -t A5SS -c {params.c} --nthread {params.thread}
        mv Darts_BHT.results.xlsx {output}
        """

rule bayes_infer_A3SS:
    input:
        rmats_count_complete="{sample}_darts/rmats_count_complete.txt",
        rmats_count_output="JC.raw.input.A3SS.txt"
    output:
        "A3SS_BHT.results.xlsx"
    params:
        annot="fromGTF.A3SS.txt",
        thread=8,
        c=0.05
    shell:
        """
        Darts_BHT bayes_infer --rmats-count {input.rmats_count_output} --od ./ --annot {params.annot} -t A3SS -c {params.c} --nthread {params.thread}
        mv Darts_BHT.results.xlsx {output}
        """

rule bayes_infer_SE:
    input:
        rmats_count_complete="{sample}_darts/rmats_count_complete.txt",
        rmats_count_output="JC.raw.input.SE.txt"
    output:
        "SE_BHT.results.xlsx"
    params:
        annot="fromGTF.SE.txt",
        thread=8,
        c=0.05
    shell:
        """
        Darts_BHT bayes_infer --rmats-count {input.rmats_count_output} --od ./ --annot {params.annot} -t SE -c {params.c} --nthread {params.thread}
        mv Darts_BHT.results.xlsx {output}
        """

rule bayes_infer_RI:
    input:
        rmats_count_complete="{sample}_darts/rmats_count_complete.txt",
        rmats_count_output="JC.raw.input.RI.txt"
    output:
        "RI_BHT.results.xlsx"
    params:
        annot="fromGTF.RI.txt",
        thread=8,
        c=0.05
    shell:
        """
        Darts_BHT bayes_infer --rmats-count {input.rmats_count_output} --od ./ --annot {params.annot} -t RI -c {params.c} --nthread {params.thread}
        mv Darts_BHT.results.xlsx {output}
        """

rule all:
    input:
        expand("alignment/{sample}.sortedByCoord.out.bam", sample=SAMPLES),
        expand("A5SS_BHT.results.xlsx", sample=SAMPLES),
        expand("A3SS_BHT.results.xlsx", sample=SAMPLES),
        expand("SE_BHT.results.xlsx", sample=SAMPLES),
        expand("RI_BHT.results.xlsx", sample=SAMPLES)

