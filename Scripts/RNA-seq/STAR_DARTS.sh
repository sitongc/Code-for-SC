#STAR to do the alignment and create the sortedbam file.
STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --runThreadN 8 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 --alignIntronMax 299999 --genomeDir ./ --sjdbGTFfile gencode.v42lift37.annotation.gtf --outFileNamePrefix ./ --readFilesIn 
       

##DART code 
Darts_BHT rmats_count --b1 b1.txt --b2 b2.txt --gtf gencode.v42lift37.annotation.gtf --od TUB0805_darts -t paired --nthread 24 --readLength 101

Darts_BHT bayes_infer --rmats-count JC.raw.input.A5SS.txt --od ./ --annot fromGTF.A5SS.txt -t A5SS -c 0.05 --nthread 8
mv Darts_BHT.results.xlsx A5SS_BHT.results.xlsx

Darts_BHT bayes_infer --rmats-count JC.raw.input.A3SS.txt --od ./ --annot fromGTF.A3SS.txt -t A3SS -c 0.05 --nthread 8
mv Darts_BHT.results.xlsx A3SS_BHT.results.xlsx

Darts_BHT bayes_infer --rmats-count JC.raw.input.SE.txt --od ./ --annot fromGTF.SE.txt -t SE -c 0.05 --nthread 8
mv Darts_BHT.results.xlsx SE_BHT.results.xlsx

Darts_BHT bayes_infer --rmats-count JC.raw.input.RI.txt --od ./ --annot fromGTF.RI.txt -t RI -c 0.05 --nthread 8
mv Darts_BHT.results.xlsx RI_BHT.results.xlsx
