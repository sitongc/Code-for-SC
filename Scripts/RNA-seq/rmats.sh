i=$(awk '{print " " $1}' data.txt)
for data in $i

do
	echo  ${data}_L001_R1_001.fastq:${data}_L001_R2_001.fastq > s1.txt
	rmats.py --s1 s1.txt --s2 s2.txt --gtf ref/Homo_sapiens.GRCh37.75.gtf --bi ./ -t paired --nthread 4 --od ${data}/ --tmp tmp_${data}/ --readLength 114 --variable-read-length 
	aws s3 cp ${data}/ s3://zuchnerlab/RNASeq/WGS_RNASeq_2022.02.24/rmats/${data}/
	sed '$d' s1.txt > s1.txt
done
