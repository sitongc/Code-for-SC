#STAR to do the alignment and create the sortedbam file.

STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --runThreadN 8 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 --alignIntronMax 299999 --genomeDir {input.index}  --sjdbGTFfile gencode.v42lift37.annotation.gtf --outFileNamePrefix ./ --readFilesIn  {input.R1} {input.R2}
        mv alignment/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}

#
