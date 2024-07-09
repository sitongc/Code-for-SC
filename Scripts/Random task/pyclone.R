library(data.table)
a <- fread('desktop/s2s39.GATK.snp.tsv')
library(tidyr)
ab<- separate(a,s2s39,into = c('ab','abc'),sep = ':')
ab<- separate(ab,abc,into = c('ref_counts','var_counts'),sep = ',')
ab$ref_counts<- as.numeric(ab$ref_counts)
ab$var_counts<- as.numeric(ab$var_counts)
ab$ddd <- ab$ref_counts +ab$var_counts
ab$AF <- ab$var_counts / ab$ddd
ab <- subset(ab,AF!=0.5)
ab <- subset(ab,AF!=0)
ab<- subset(ab,AF!=1)
ab$mutation_id <-paste(ab$`#CHROM`,ab$POS,sep = ':')
ab$normal_cn <- 2
ab$minor_cn <- 0
ab$major_cn <- 2
abb <- ab[,11:18]
abb<- abb[,-3:-4]
abb <- abb[,c(3,1,2,4,5,6)]
abb <- na.omit(abb)
write.table(abb,'desktop/s2s39.tsv',quote=FALSE, sep='\t',row.names = FALSE)
