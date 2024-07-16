#!/bin/bash
# jointVariantCallGVCFs_v1

## 
# INPUT: cohortRunID -- the GUID assigned to this run of the cohort pipeline
# Query metastore to get the number of samples and the location of their gVCF files
# Perform JVC on the gVCF files
# Change the status of this family to joint-called and upload the output file to S3

set -euo pipefail

source aws-pipeline/utils/utils.sh  # imports errorReport

sqlcmd -S tgp-metastore.cqbat9wvtmzz.us-east-1.rds.amazonaws.com,1433 -U development -P TGPdeveloper! -d prod -Q "select @@version" > /dev/null 2>&1 || true

cohortRunID=$1
cohortRunID=${cohortRunID,,} # ensure that the ID is lowercase

trap "errorReport cohort ${cohortRunID}" ERR

cd /data

# input parameters
gatk="/data/gatk-4.1.4.1/gatk"
intervals="/data/aws-pipeline/jointVariantCalling/allCallingIntervals.list"
ref_fasta="/data/hs37d5.fa"
dbsnp_vcf="/data/ref/dbsnp_138.b37.vcf"

cp aws-pipeline/jointVariantCalling/getSampleLocations.sql .
sed -i "s/<cohortRunID>/${cohortRunID}/" getSampleLocations.sql


# get the s3 locations of each gVCF file and the file name
# this will return a text file table called metastoreResponse.txt 
sqlcmd -S tgp-metastore.cqbat9wvtmzz.us-east-1.rds.amazonaws.com,1433 -U development -P TGPdeveloper! -d prod -i getSampleLocations.sql -o metastoreResponse.txt
# getSampleLocations.sql selects the the sampleID, sample-level run ID, gVCF file name, and upload type of each member of the queried family
# groom the metastoreResponse.txt file a little bit
tail -n +3 metastoreResponse.txt | head -n -2 > tmp; mv tmp metastoreResponse.txt 

#update the pipeline to say that this cohort has started 
sqlcmd -S tgp-metastore.cqbat9wvtmzz.us-east-1.rds.amazonaws.com,1433 -U development -P TGPdeveloper! -d prod -Q "exec UpdateCohortRun @CohortRunId = '${cohortRunID}', @Status = 'Running'"

# transfer the gVCF files to this workspace
GVCFFiles=""
#anyGenomes="false"
while read sampleID runID fileName uploadType sampleGUID; do 
	fileName=$(sed 's/\.gz//' <<< ${fileName,,})
	aws s3 cp ${fileName}.gz . --only-show-errors
	aws s3 cp ${fileName}.idx . --only-show-errors
	BASE=$(sed 's#^.*/##g' <<< ${fileName})
	if [ ! -f ${BASE} ]; then 
		gunzip ${BASE}.gz
	fi
	GVCFFiles+=" -V "${BASE}
done < metastoreResponse.txt 

# clean up the GVCFFiles list
GVCFFiles=$(sed 's/^-V //' <<< ${GVCFFiles}) # personal preference
GVCFFiles=$(sed 's/ -V$//' <<< ${GVCFFiles}) # just in case there is a trailing empty line

# calculate how much memory is available (this will round down which underestimates a bit, but that is desired in this case)
MEMKB=$(cat /proc/meminfo | grep 'MemTotal' | sed 's/[ ]\+/\t/g' | cut -f 2)
MEMMB=$(( MEMKB / 1024 ))
MEMGB=$(( MEMMB / 1024 ))

# prepare 1000Genomes data on allele frequencies and genetic mappings for ROH mapping
aws s3 cp s3://tgp-data-analysis/StructuralVariation/genetic-map.tgz . --only-show-errors
aws s3 cp s3://tgp-data-analysis/StructuralVariation/1000GP-AFs.tgz . --only-show-errors
tar -zxvf 1000GP-AFs.tgz
tar -zxvf genetic-map.tgz

# Perform CNN filtering on each individually-called VCF
set +eu
while read sampleID runID fileName uploadType sampleGUID; do 
	fileName=$(sed 's/\.gz//' <<< ${fileName,,})
	fileName=$(sed 's#^.*/##g' <<< ${fileName})
	BASE=$(sed 's/\.hg19\.g\.vcf//' <<< ${fileName})
	${gatk} --java-options "-Xmx${MEMGB}g" GenotypeGVCFs -V ${fileName} -R ${ref_fasta} -O ${BASE}.hg19.genotyped.vcf --allow-old-rms-mapping-quality-annotation-data --use-new-qual-calculator --only-output-calls-starting-in-intervals -L ${intervals}
	bcftools filter -i "(FORMAT/DP)>=6" ${BASE}.hg19.genotyped.vcf > ${BASE}.hg19.genotyped.filtered.vcf
	source activate gatk
	${gatk} CNNScoreVariants -V ${BASE}.hg19.genotyped.filtered.vcf -R ${ref_fasta} -O ${BASE}.CNN.vcf
	${gatk} FilterVariantTranches -V ${BASE}.CNN.vcf -O ${BASE}.CNN.filtered.vcf --resource /data/ref/hapmap_3.3.b37.vcf --resource /data/ref/Mills_and_1000G_gold_standard.indels.b37.vcf --info-key CNN_1D --snp-tranche 99.95 --snp-tranche 99.9 --snp-tranche 99.5 --snp-tranche 99 --snp-tranche 98 --snp-tranche 95 --indel-tranche 99.4 --indel-tranche 99.2 --indel-tranche 99 --indel-tranche 98 --indel-tranche 95
	conda deactivate
	# upload the CNN-filtered VCF to the cohortRun folder on S3
	aws s3 cp ${BASE}.CNN.filtered.vcf s3://tgp-sample-processing/${cohortRunID}/ --only-show-errors
	## Runs of Homozygosity mapping
	# remove the CNN-failed variants
	grep -wv 'CNN_1D_SNP_Tranche_99.95_100.00' ${BASE}.CNN.filtered.vcf | grep -wv 'CNN_1D_INDEL_Tranche_99.40_100.00' > ${BASE}.filtered.vcf
	bgzip ${BASE}.filtered.vcf
	tabix ${BASE}.filtered.vcf.gz
	# run bcftools roh on the sample
	bcftools annotate -c CHROM,POS,REF,ALT,AF1KG -h 1000GP-AFs/AFs.tab.gz.hdr -a 1000GP-AFs/AFs.tab.gz ${BASE}.filtered.vcf.gz | bcftools roh --AF-tag AF1KG -M 100 -m genetic-map/genetic_map_chr{CHROM}_combined_b37.txt -o roh.txt
	grep '^RG' roh.txt > ${BASE}.RohRegions.txt
	# filter the output to only include ROHs at least 5kb in size
	echo -e "chr\tstart\tend\tlength\tnumberOfMarkers\tquality" > header.txt
	awk -F'\t' -v OFS='\t' '{if ($6>5000){print $3,$4,$5,$6,$7,$8}}' ${BASE}.RohRegions.txt | cat header.txt - > ${BASE}.RegionsOfHomozygosity.txt
	# create bed file for IGV visualization
	echo "#gffTags" > header.txt
	awk -F'\t' -v OFS='\t' '{if ($6>5000){print $3,$4,$5,"Length="$6";NumberOfMarkers="$7";Quality="$8}}' ${BASE}.RohRegions.txt | cat header.txt - > ${BASE}.RegionsOfHomozygosity.bed
	# upload the roh's back to s3
	aws s3 cp ${BASE}.RegionsOfHomozygosity.txt s3://tgp-sample-processing/${runID,,}/ --only-show-errors
	aws s3 cp ${BASE}.RegionsOfHomozygosity.bed s3://tgp-sample-assets/${sampleGUID,,}/ROHs/ --only-show-errors
	rm ${BASE}.hg19.genotyped.vcf ${BASE}.CNN.vcf* ${BASE}.CNN.filtered.vcf* ${BASE}.filtered* ${BASE}.RohRegions.txt ${BASE}.RegionsOfHomozygosity* roh.txt
done < metastoreResponse.txt 
set -eu


# Combine the gVCFs using combineGVCFs
${gatk} --java-options "-Xmx${MEMGB}g" CombineGVCFs \
	-R ${ref_fasta} \
	-V ${GVCFFiles} \
	-O ${cohortRunID}.combined.g.vcf


# joint genotype the gVCFs
${gatk} --java-options "-Xmx${MEMGB}g" GenotypeGVCFs \
	-R ${ref_fasta} \
	-O ${cohortRunID}.vcf \
	--allow-old-rms-mapping-quality-annotation-data \
	--only-output-calls-starting-in-intervals \
	--use-new-qual-calculator \
	-V ${cohortRunID}.combined.g.vcf \
	-L ${intervals}


# filter the output file
FILE="${cohortRunID}.vcf"
chmod 777 aws-pipeline/jointVariantCalling/bcftools
# include only sites covered in at least one sample
aws-pipeline/jointVariantCalling/bcftools filter -i "MAX(FORMAT/DP)>=6" ${FILE} > tmp
mv tmp ${FILE}
# include only sites that do not have anomalously high coverage in all samples (unless all samples are panels) --removed
#if [ "${allPanels}" == "false" ]; then
#	aws-pipeline/jointVariantCalling/bcftools filter -i "MIN(FORMAT/DP)<500" tmp > ${FILE}
#else
#	mv tmp ${FILE}
#fi

# filter vcf file to only include the exome capture regions used by ExAC -- removed 
#if [ "${anyGenomes}" == "false" ]; then
#	./bedtools intersect -u -header -a ${FILE} -b exome_calling_regions.bed > tmp
#	mv tmp ${FILE}
#fi

# perform sample-wise vcfLevelQC
perlPath="/data/perl/"
vcfTSTV="/data/vcf-tstv"
vcfStats="/data/vcf-stats"
vcftools="/data/vcftools"
vcfSubset="/data/perl/vcf-subset"
export PERL5LIB=${perlPath}

while read sampleID runID fileName uploadType sampleGUID; do
	BASE=${sampleID}
	${vcfSubset} -c ${BASE} ${FILE} > ${BASE}.vcf
	echo -e "fileName\tsampleName\tChr\tTSTV\tnumSNPs\tnumINDELs\tavgDepth\tavgQuality\tHeterozygousHomozygousRatio\tmeanGQ\tGender" > ${BASE}.vcfLevel.byChr.qc
	TSTV=$(cat ${BASE}.vcf | ${vcfTSTV} | awk '{print $1}')
	${vcfStats} ${BASE}.vcf > tmpStats.txt 
	SNPs=$((grep 'snp_count' tmpStats.txt || :) | sed '1q;d' | awk '{print $3}' | sed 's/,$//') 
	INDELs=$((grep 'indel_count' tmpStats.txt || :) | sed '1q;d' | awk '{print $3}' | sed 's/,$//') 
	DEPTH=$(${vcftools} --vcf ${BASE}.vcf --depth --stdout | sed '2q;d' | awk '{print $3}')
	QUAL=$(${vcftools} --vcf ${BASE}.vcf --site-quality --stdout | tail -n +2 | awk '{total += $3} END {if (NR>0) {print total/NR} else {print "NA"}}')
	HET_RA=$((grep 'het_RA' tmpStats.txt || :) | sed '1q;d' | awk '{print $3}' | sed 's/,$//')
	HET_AA=$((grep 'het_AA' tmpStats.txt || :) | sed '1q;d' | awk '{print $3}' | sed 's/,$//')
	HOM_AA=$((grep 'hom_AA' tmpStats.txt || :) | sed '1q;d' | awk '{print $3}' | sed 's/,$//')
	if [ "${SNPs}" == "" ]; then SNPs="0"; fi
	if [ "${INDELs}" == "" ]; then INDELs="0"; fi
	if [ "${HET_RA}" == "" ]; then HET_RA="0"; fi
	if [ "${HET_AA}" == "" ]; then HET_AA="0"; fi
	if [ "${HOM_AA}" == "" ]; then HOM_AA="0"; fi
	HETHOMRATIO=$(awk -v het1="${HET_RA}" -v het2="${HET_AA}" -v hom="${HOM_AA}" 'BEGIN {if (hom>0) {print (het1+het2)/hom} else {print "NA"}}')
	GQ=$(${vcftools} --vcf ${BASE}.vcf --extract-FORMAT-info GQ --stdout | tail -n +2 | awk '{sum+=$3} END {if (NR>0) {print sum/NR} else {print "NA"}}')
	## old, plink-based sex check. Kept for reference
	#./plink --vcf ${BASE}.vcf --split-x b37 no-fail --make-bed --out ${BASE}
	#./plink --bfile ${BASE} --check-sex --out ${BASE}
	#SEX=$(sed '2q;d' ${BASE}.sexcheck | awk '{print $4}')
	#if [ ${SEX} == 1 ]; then
	#	SEX='M'
	#elif [ ${SEX} == 2 ]; then
	#	SEX='F'
	#else 
	#	SEX='U'
	#fi 
	## new, peddy-based sex check
	echo -e "${BASE}\t${BASE}\t0\t0\t0\t0" > ${BASE}.ped
	bgzip -c ${BASE}.vcf > ${BASE}.vcf.gz
	tabix ${BASE}.vcf.gz
	/home/ec2-user/anaconda/bin/python -m peddy --prefix ${BASE} ${BASE}.vcf.gz ${BASE}.ped
	SEX=$(sed -n '2p' ${BASE}.sex_check.csv | cut -f 7 -d, )
	if [ ${SEX} == "male" ]; then
		SEX='M'
	elif [ ${SEX} == "female" ]; then
		SEX='F'
	else 
		SEX='U'
	fi 
	
	echo -e "${FILE}\t${BASE}\tGenome\t${TSTV}\t${SNPs}\t${INDELs}\t${DEPTH}\t${QUAL}\t${HETHOMRATIO}\t${GQ}\t${SEX}" >> ${BASE}.vcfLevel.byChr.qc 
	for j in {1..22} X Y; do 
		${vcftools} --vcf ${BASE}.vcf --chr ${j} --stdout --recode > thisChr.vcf 
		numVariants=$(grep -v ^'#' thisChr.vcf | wc -l)
		if [ "${numVariants}" -gt "0" ]; then 
			TSTV=$(cat thisChr.vcf | ${vcfTSTV} | awk '{print $1}')
			${vcfStats} thisChr.vcf > tmpStats.txt 
			SNPs=$((grep 'snp_count' tmpStats.txt || :) | sed '1q;d' | awk '{print $3}' | sed 's/,$//') 
			INDELs=$((grep 'indel_count' tmpStats.txt || :) | sed '1q;d' | awk '{print $3}' | sed 's/,$//') 
			DEPTH=$(${vcftools} --vcf thisChr.vcf --depth --stdout | sed '2q;d' | awk '{print $3}')
			QUAL=$(${vcftools} --vcf thisChr.vcf --site-quality --stdout | tail -n +2 | awk '{total += $3} END {if (NR>0) {print total/NR} else {print "NA"}}')
			HET_RA=$((grep 'het_RA' tmpStats.txt || :) | sed '1q;d' | awk '{print $3}' | sed 's/,$//')
			HET_AA=$((grep 'het_AA' tmpStats.txt || :) | sed '1q;d' | awk '{print $3}' | sed 's/,$//')
			HOM_AA=$((grep 'hom_AA' tmpStats.txt || :) | sed '1q;d' | awk '{print $3}' | sed 's/,$//')
			if [ "${SNPs}" == "" ]; then SNPs="0"; fi
			if [ "${INDELs}" == "" ]; then INDELs="0"; fi
			if [ "${HET_RA}" == "" ]; then HET_RA="0"; fi
			if [ "${HET_AA}" == "" ]; then HET_AA="0"; fi
			if [ "${HOM_AA}" == "" ]; then HOM_AA="0"; fi
			HETHOMRATIO=$(awk -v het1="${HET_RA}" -v het2="${HET_AA}" -v hom="${HOM_AA}" 'BEGIN {if (hom>0) {print (het1+het2)/hom} else {print "NA"}}')
			GQ=$(${vcftools} --vcf thisChr.vcf --extract-FORMAT-info GQ --stdout | tail -n +2 | awk '{sum+=$3} END {if (NR>0) {print sum/NR} else {print "NA"}}')
			echo -e "${FILE}\t${BASE}\t${j}\t${TSTV}\t${SNPs}\t${INDELs}\t${DEPTH}\t${QUAL}\t${HETHOMRATIO}\t${GQ}\t${SEX}" >> ${BASE}.vcfLevel.byChr.qc 
		fi
	done 
	##ancestry clustering
	runID=${runID,,}
	aws s3 cp s3://tgp-sample-processing/${runID}/${BASE}.peddySites.vcf . --only-show-errors
	${vcftools} --vcf ${BASE}.peddySites.vcf --012 --out ${BASE}
	/usr/bin/python samplePCA.py -i ${BASE}
	#upload results 
	aws s3 cp ${BASE}.AncestryClustering.tiff s3://tgp-sample-processing/${cohortRunID}/ --only-show-errors
	aws s3 cp ${BASE}.vcfLevel.byChr.qc s3://tgp-sample-processing/${cohortRunID}/ --only-show-errors
	# tag the cohortRun output files so they can be archived and deleted as appropriate
	EXISTS=$(aws s3 ls s3://tgp-sample-processing/${cohortRunID}/${BASE}.vcfLevel.byChr.qc | grep ${BASE} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3api put-object-tagging --bucket tgp-sample-processing --key ${cohortRunID}/${BASE}.vcfLevel.byChr.qc --tagging '{ "TagSet": [{ "Key": "Archive", "Value": "True"}]}'; fi
	EXISTS=$(aws s3 ls s3://tgp-sample-processing/${cohortRunID}/${BASE}.AncestryClustering.tiff | grep ${BASE} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3api put-object-tagging --bucket tgp-sample-processing --key ${cohortRunID}/${BASE}.AncestryClustering.tiff --tagging '{ "TagSet": [{ "Key": "Archive", "Value": "True"}]}'; fi
	# tag the sampleRun output files so they can be archived and deleted as appropriate
	EXISTS=$(aws s3 ls s3://tgp-sample-processing/${runID}/${BASE}.hg19.g.vcf.gz | grep ${BASE} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3api put-object-tagging --bucket tgp-sample-processing --key ${runID}/${BASE}.hg19.g.vcf.gz --tagging '{ "TagSet": [{ "Key": "Archive", "Value": "True"}]}'; fi
	EXISTS=$(aws s3 ls s3://tgp-sample-processing/${runID}/${BASE}.hg19.g.vcf.idx | grep ${BASE} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3api put-object-tagging --bucket tgp-sample-processing --key ${runID}/${BASE}.hg19.g.vcf.idx --tagging '{ "TagSet": [{ "Key": "Archive", "Value": "True"}]}'; fi
	EXISTS=$(aws s3 ls s3://tgp-sample-processing/${runID}/${BASE}.peddySites.vcf | grep ${BASE} | wc -l || true)
	#if [[ ${EXISTS} -ne 0 ]]; then aws s3api put-object-tagging --bucket tgp-sample-processing --key ${runID}/${BASE}.peddySites.vcf --tagging '{ "TagSet": [{ "Key": "Archive", "Value": "True"}]}'; fi
	EXISTS=$(aws s3 ls s3://tgp-sample-processing/${runID}/${BASE}.hg19.bamLevel.qc | grep ${BASE} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3api put-object-tagging --bucket tgp-sample-processing --key ${runID}/${BASE}.hg19.bamLevel.qc --tagging '{ "TagSet": [{ "Key": "Archive", "Value": "True"}]}'; fi
	EXISTS=$(aws s3 ls s3://tgp-sample-processing/${runID}/${BASE}.fastqLevel.qc | grep ${BASE} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3api put-object-tagging --bucket tgp-sample-processing --key ${runID}/${BASE}.fastqLevel.qc --tagging '{ "TagSet": [{ "Key": "Archive", "Value": "True"}]}'; fi
	EXISTS=$(aws s3 ls s3://tgp-sample-processing/${runID}/${BASE}.sorted.bam | grep ${BASE} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3api put-object-tagging --bucket tgp-sample-processing --key ${runID}/${BASE}.sorted.bam --tagging '{ "TagSet": [{ "Key": "Archive", "Value": "False"}]}'; fi
	EXISTS=$(aws s3 ls s3://tgp-sample-processing/${runID}/${BASE}.hg19.sorted.dedup.bam | grep ${BASE} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3api put-object-tagging --bucket tgp-sample-processing --key ${runID}/${BASE}.hg19.sorted.dedup.bam --tagging '{ "TagSet": [{ "Key": "Archive", "Value": "False"}]}'; fi
	EXISTS=$(aws s3 ls s3://tgp-sample-processing/${runID}/${BASE}.hg19.sorted.dedup.bam.bai | grep ${BASE} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3api put-object-tagging --bucket tgp-sample-processing --key ${runID}/${BASE}.hg19.sorted.dedup.bam.bai --tagging '{ "TagSet": [{ "Key": "Archive", "Value": "False"}]}'; fi
	EXISTS=$(aws s3 ls s3://tgp-sample-processing/${runID}/${BASE}.RegionsOfHomozygosity.txt | grep ${BASE} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3api put-object-tagging --bucket tgp-sample-processing --key ${runID}/${BASE}.RegionsOfHomozygosity.txt --tagging '{ "TagSet": [{ "Key": "Archive", "Value": "True"}]}'; fi
done < metastoreResponse.txt 

# upload the results back to S3
gzip ${cohortRunID}.vcf
aws s3 cp ${cohortRunID}.vcf.gz s3://tgp-sample-processing/${cohortRunID}/ --only-show-errors

# transfer the input fastq files to the archiving S3 account
cp aws-pipeline/jointVariantCalling/getInputFileLocations.sql .
sed -i "s/<cohortRunID>/${cohortRunID}/" getInputFileLocations.sql
sqlcmd -S tgp-metastore.cqbat9wvtmzz.us-east-1.rds.amazonaws.com,1433 -U development -P TGPdeveloper! -d prod -i getInputFileLocations.sql -o inputFileLocations.txt
while read FileLink1 FileLink2; do 
	FileLink1=$(sed 's#s3://tgp-s3-bucket-ftp/##' <<< ${FileLink1})
	FileLink2=$(sed 's#s3://tgp-s3-bucket-ftp/##' <<< ${FileLink2})
	FileName1=$(sed 's#^.*/##' <<< ${FileLink1})
	FileName2=$(sed 's#^.*/##' <<< ${FileLink2})
	EXISTS=$(aws s3 ls "s3://tgp-s3-bucket-ftp/${FileLink1}" | grep ${FileName1} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3 cp "s3://tgp-s3-bucket-ftp/${FileLink1}" "s3://genesis-files-to-archive/${FileLink1}" --acl bucket-owner-full-control --only-show-errors; aws s3 rm s3://tgp-s3-bucket-ftp/${FileLink1}; fi
	EXISTS=$(aws s3 ls "s3://tgp-s3-bucket-ftp/${FileLink2}" | grep ${FileName2} | wc -l || true)
	if [[ ${EXISTS} -ne 0 ]]; then aws s3 cp "s3://tgp-s3-bucket-ftp/${FileLink2}" "s3://genesis-files-to-archive/${FileLink2}" --acl bucket-owner-full-control --only-show-errors; aws s3 rm s3://tgp-s3-bucket-ftp/${FileLink2}; fi
done < inputFileLocations.txt


#update the pipeline to say that this sample has finished -- no longer needed with addition of parser module
#sqlcmd -S tgp-metastore.cqbat9wvtmzz.us-east-1.rds.amazonaws.com,1433 -U development -P TGPdeveloper! -d prod -Q "exec UpdateCohortRun @CohortRunId = '${cohortRunID}', @Status = 'Done'"

