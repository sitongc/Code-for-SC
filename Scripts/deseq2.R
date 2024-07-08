dir <- "C:/Users/Sitong Chen/Desktop/salmon/"
list.files(dir)
samples <- paste0("Sample", c(seq(1,6)))
files <- file.path(dir, samples, "quant.sf")
##sample name
names(files) <-paste0(c('iCN SP-63 3.1 A','iCN SP-64 3.1 A','iCN AX-45 isog. Control #25 A','F-AX-45 B','F-SP-64','F-SP-63'))
all(file.exists(files))

# reference gene-level annotation package or input the reference file.
library(EnsDb.Hsapiens.v86)
txdb <- EnsDb.Hsapiens.v86
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# associated salmon output with annotation package
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE) 
names(txi)
txi$length[txi$length == 0] <- 1

#Using DESeq2 to analysis differential gene expression
library(DESeq2)
sampleTable <- data.frame(condition = factor(rep(c("control", "patient"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, name="condition_patient_vs_control")
res <- results(dds, contrast=c("condition","patient","control"))
res

#Log fold change shrinkage(reduce noises)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_patient_vs_control", type="apeglm")
resLFC

#p-values and adjusted p-values
resOrdered <- res[order(res$pvalue),] #order our results table by the smallest p value
summary(res)
sum(res$padj < 0.1, na.rm=TRUE) #how many adjusted p-values were less than 0.1
res01 <- results(dds, alpha=0.1) #default alpha is set to 0.1, if the adjusted p value cutoff will be a value other than 0.1, alpha = 0.05
summary(res01)

#MA-plot
plotMA(res, ylim=c(-2,2))

#Plot counts: examine the counts of reads for a single gene across the group 
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#Export results to CVS files
write.csv(as.data.frame(resOrdered), 
          file="HSP.csv")

#Extracting transformed values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

##Extracting transformed values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

#Effects of transformations on the variance (standard deviation of transformed data, across samples, aginst the mean, using the shifted logarithm transformation)library("vsn"))
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

#heatmap: various transformation of data, choose the specific gene you are interest in. 
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20] #stand for degree of color
df <- as.data.frame(colData(dds)[c("condition")])
rownames(df) <- colnames(dds)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)#top 20 gene.

#heatmap of the sample to sample distances
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

#principal component plot of the samples
plotPCA(vsd, intgroup=c("condition"))

#count outliers
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#dispersion plot
plotDispEsts(dds)