# Code-for-Sitong Chen


## RNA-seq pipeline
![example output](RNA-pipeline.png)

All of the scripts are available in **Scripts/RNA-seq/**. 

Besides the softwares were mentions in the pipeline, this part also include 2 other splicing variant softwares, rMATS and Leafcutter. 

It also provided the **Snakemake** version. 





## Linux 

All of the process running in HPC

### Variance calling
#### BWA mapping 
```
for id in /groups/ligrp/sitongchen/scrna/sc/*_1.fastq.gz
do
date
filename=$(basename ${id} _1.fastq.gz)
filepath=${id%/*}; echo "BWA mapping for ${filename}"
bwa mem -M -t 20 /groups/ligrp/sitongchen/scrna/22/GRCm38_68.fa ${filepath}/${filename}_1.fastq.gz ${filepath}/${filename}_2.fastq.gz > groups/ligrp/sitongchen/scrna/dpi35/result/sam/${filename}.sam
```
#### Sort sam files

```
date
echo "Sort sam file"
gatk SortSam \
--INPUT /groups/ligrp/sitongchen/scrna/sam/${filename}.sam \
--OUTPUT /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted.bam \
--SORT_ORDER coordinate
if [ /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted.bam -s ]; then
rm -f /groups/ligrp/sitongchen/scrna/sam/${filename}.sam
fi
```

#### Marking duplicates
```
date
echo "Marking duplicates"
gatk MarkDuplicates \
--INPUT /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted.bam \
--METRICS_FILE /groups/ligrp/sitongchen/scrna/sortbam/${filename}_dupMetrics \
--OUTPUT /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted_dupMarked.bam
if [ /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted_dupMarked.bam -s ]; then
rm -f /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted.bam
fi
```
##### Repairing Readgroups
```
    date
    echo "Repairing Readgroups";
    gatk AddOrReplaceReadGroups \
    --INPUT /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted_dupMarked.bam \
    --OUTPUT /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted_DM_RG.bam \
    --RGLB illumina_${filename} \
    --RGPL illumina \
    --RGPU JR001 \
    --RGSM ${filename}
```

#### First recalibration table

```
gatk BaseRecalibrator \
--input /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted_DM_RG.bam \
--known-sites /groups/ligrp/sitongchen/scrna/22/mgp.v3.snps.rsIDdbSNPv137.vcf.gz \
--output /groups/ligrp/sitongchen/scrna/sortbam/${filename}_recal.table \
--reference /groups/ligrp/sitongchen/scrna/22/GRCm38_68.fa

gatk ApplyBQSR \
        -R /groups/ligrp/sitongchen/scrna/22/GRCm38_68.fa \
        -I  /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted_DM_RG.bam \
        --bqsr-recal-file /groups/ligrp/sitongchen/scrna/sortbam/${filename}_recal.table \
        -O  /groups/ligrp/sitongchen/scrna/sortbam/${filename}_BR.bam
```
#### VCF file generation
```
gatk Mutect2 \
-R /groups/ligrp/sitongchen/scrna/22/GRCm38_68.fa \
-I /groups/ligrp/sitongchen/scrna/sortbam/${filename}_BR.bam \
-O /groups/ligrp/sitongchen/scrna/gatk/${filename}_raw.vcf.gz

```

#### BCFtool Variance calling 
```
for id in /groups/ligrp/sitongchen/scrna/dpi35/result/sortbam/*_BR.bam
do
    filename=$(basename ${id} _BR.bam)
    echo "  -->  process:  $filename"
    bcftools mpileup --threads 10 -q 20 -Ou -f /groups/ligrp/sitongchen/scrna/22/GRCm38_68.fa  /groups/ligrp/sitongchen/scrna/dpi35/result/sortbam/${filename}_BR.bam | bcftools call --threads 4 -mv -Ov  > /groups/ligrp/sitongchen/scrna/dpi35/result/bcfs/${filename}.raw.vcf
done
```

#### scSplitter
scSplitter is a preprocessing tool designed to convert scRNA seq results into a format suitable for mutation or SNV calling at the single cell level for the purpose of lineage tracing.

```
STAR --runMode genomeGenerate --genomeDir /groups/ligrp/sitongchen/scrna/dpi35  --genomeFastaFiles /groups/ligrp/sitongchen/scrna/22/mm10.fa --sjdbGTFfile /groups/ligrp/sitongchen/scrna/22/mm10.refGene.gtf -- sjdbOverhang 199

python3 scSplitter.py --ig True --f /groups/ligrp/sitongchen/scrna/dpi35/gel_barcode1_list.txt --i /groups/ligrp/sitongchen/scrna/dpi35/fastqs --r /groups/ligrp/sitongchen/scrna/dpi35/result --ind /groups/ligrp/sitongchen/scrna/dpi35/sta
STAR --runMode genomeGenerate --genomeDir /groups/ligrp/sitongchen/scrna/dpi35  --genomeFastaFiles /groups/ligrp/sitongchen/scrna/22/mm10.fa --sjdbGTFfile /groups/ligrp/sitongchen/scrna/22/mm10.refGene.gtf -- sjdbOverhang 199

```

#### GATK liftOver

liftOver a VCF file from one reference build to another.

In my case, I will liftover GRCh38 to hg19. You need to download the hg19.fa and hg38ToHg19.over.chain file for this step. Based on the size of the ref, you may need a instance with large RAM. 

```
java -Xmx370G -jar picard.jar  LiftoverVcf --CHAIN hg38ToHg19.over.chain --INPUT gnomad.genomes.v4.1.sites.chr${data}.vcf  --OUTPUT chr${data}.vcf --REFERENCE_SEQUENCE hg19.fa --REJECT reject.vcf
```


