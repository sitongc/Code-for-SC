i=$(awk '{print " " $1}' sample.txt)
for BASE in $i
do
aws s3 cp s3://zuchnerlab/testCases/ . --recursive --exclude "*" --include "${BASE}*" 
# change the name of the bam and bai files to fit the expected pattern
EXISTS=$(ls ${BASE}.hg19.sorted.dedup.bam | grep ${BASE} | wc -l || true)
if [[ ${EXISTS} -eq 0 ]]; then
	mv ${BASE}*.bam ${BASE}.hg19.sorted.dedup.bam
	mv ${BASE}*.bai ${BASE}.hg19.sorted.dedup.bam.bai
fi

echo "`date`"

ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter --reads ${BASE}.hg19.sorted.dedup.bam --reference ref/hs37d5.fa --variant-catalog variant_catalog.json --output-prefix ${BASE}
samtools sort ${BASE}_realigned.bam > ${BASE}_realigned.sorted.bam
samtools index ${BASE}_realigned.sorted.bam
grep -v '^#' ${BASE}.vcf | cut -f 5 -d';' | sed 's/VARID=//' | grep -v '_' > varIDs.txt

while read VARID; do
	./REViewer-v0.2.7-linux_x86_64 --reads ${BASE}_realigned.sorted.bam --vcf ${BASE}.vcf --reference ref/hs37d5.fa --catalog variant_catalog.json --locus ${VARID} --output-prefix ${BASE}
done < varIDs.txt

echo "`date`"

done

