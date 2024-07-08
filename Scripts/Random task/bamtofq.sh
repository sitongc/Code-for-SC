##samtools 1.12

##samtools collate

samtools collate -o name-collate.bam 3cd1e613-829b-40dd-aa6a-3460085eec07_wgs_gdc_realn.bam
samtools fastq -1 sample_1.fastq.gz -2 sample_2.fastq.gz  name-collate.bam


##bedtools

sambamba sort -m 5GB -t 2 -n -o ${BASE}.nameSorted.bam ${FILE}
rm ${FILE}
/data/bedtools bamtofastq -i ${BASE}.nameSorted.bam -fq ${BASE}_R1.fastq -fq2 ${BASE}_R2.fastq

##samtools bam2fq

samtools bam2fq 3cd1e613-829b-40dd-aa6a-3460085eec07_wgs_gdc_realn.bam > SAMPLE.fastq

cat SAMPLE.fastq | grep '^@.*/1$' -A 3 --no-group-separator > SAMPLE_R1.fastq

cat SAMPLE.fastq | grep '^@.*/2$' -A 3 --no-group-separator > SAMPLE_R2.fastq

##samtools sort and fastq
samtools sort  -o SAMPLE_sorted.bam 3cd1e613-829b-40dd-aa6a-3460085eec07_wgs_gdc_realn.bam
samtools fastq -1 sample_1.fastq.gz -2 sample_2.fastq.gz SAMPLE_sorted.bam