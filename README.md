# Code-for-SC
The code for R &amp; Linux &amp; Python

## R

### DESeq2

#### Import salmon results
```

dir <- "Documents/lab/DEseq2/quants"
list.files(dir)
samples <- paste0("Mouse2_", c(seq(1,4)),"_quant")
files <- file.path(dir, samples, "quant.sf")
names(files) <-paste0(c('Mouse 2-SVZL-P3','Mouse 2-SVZR-P2','R172H SVZ L4P2','R172H SVZ L5P2','R172H SVZ L6P2','181004#1 SVZ P2','181004#4 SVZ P2','181004#5 SVZ P2','181004#2 SVZ P2','181004#3 SVZ P2')) # sample name
all(file.exists(files))
```
#### Reference gene-level annotation package
```
library(EnsDb.Mmusculus.v79)
txdb <- EnsDb.Mmusculus.v79
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
```

#### Associated salmon output with annotation package

```
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
names(txi)
```

#### Using DESeq2 to analysis differential gene expression
```
library(DESeq2)
sampleTable <- data.frame(condition = factor(c(rep("SVZ",2), rep("Control", 8)))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
```

#### Differential expression analysis
```
dds <- DESeq(dds)
res <- results(dds, name="condition_treated_vs_untreated")
res <- results(dds, contrast=c("condition","treated","untreated"))
res
```

#### Log fold change shrinkage(reduce noises)
```
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC
```

#### p-values and adjusted p-values
```
resOrdered <- res[order(res$pvalue),] 
summary(res)
sum(res$padj < 0.1, na.rm=TRUE) 
res01 <- results(dds, alpha=0.1) 
summary(res01)
```

#### MA-plot
Plot counts: examine the counts of reads for a single gene across the group 
```
plotMA(res, ylim=c(-2,2))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
```

#### Export results to CVS files
```
write.csv(as.data.frame(resOrdered), 
          file="desktop/SVZ_control.csv")
```
          
#### Extracting transformed values
```
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
```

#### Effects of transformations on the variance
standard deviation of transformed data, across samples, aginst the mean, using the shifted logarithm transformation)
```
ntd <- normTransform(dds)

library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))
```

#### Heatmap: various transformation of data
Based on the DEG analysis, we can create a heatmap to visualize the data or chose the gene you interested in. Then Convert Ensembl ID to gene name. 
```
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20] # The 20 genes with the highest expression levels
df <- as.data.frame(colData(dds)[c("condition")])
rownames(df) <- colnames(dds)
pheatmap(assay(ntd)[c('ENSMUSG00000000093',
                      'ENSMUSG00000000103',
                      'ENSMUSG00000000305'),], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df) # the gene you are interested in. 

mat <- assay(ntd)[c('ENSMUSG00000000093',
                      'ENSMUSG00000000103',
                      'ENSMUSG00000000305'),]
                      
library(biomaRt)
mart <- useMart("ensembl","mmusculus_gene_ensembl", host = 'uswest.ensembl.org', ensemblRedirect = FALSE)
gns <- getBM(c("external_gene_name","ensembl_gene_id"), "ensembl_gene_id", row.names(mat), mart)
row.names(mat)[match(gns[,2], row.names(mat))] <- gns[,1] #notice the order of geneiD and gene name
pheatmap(mat, show_rownames=TRUE, annotation_col=df,display_numbers =TRUE)        
```

#### heatmap of the sample to sample distances
```
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```
#### principal component plot of the samples
```
plotPCA(vsd, intgroup=c("condition"))
```
#### count outliers
```
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
```
#### dispersion plot
```
plotDispEsts(dds)
```
