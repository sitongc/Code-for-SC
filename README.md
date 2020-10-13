# Code-for-SC
The code for R &amp; Linux &amp; Python

# R
# DESeq2
# import salmon results
dir <- "Documents/lab/DEseq2/quants"
list.files(dir)
samples <- paste0("Mouse2_", c(seq(1,4)),"_quant")
files <- file.path(dir, samples, "quant.sf")
names(files) <-paste0(c('Mouse 2-SVZL-P3','Mouse 2-SVZR-P2','R172H SVZ L4P2','R172H SVZ L5P2','R172H SVZ L6P2','181004#1 SVZ P2','181004#4 SVZ P2','181004#5 SVZ P2','181004#2 SVZ P2','181004#3 SVZ P2')) # sample name
all(file.exists(files))

# reference gene-level annotation package
library(EnsDb.Mmusculus.v79)
txdb <- EnsDb.Mmusculus.v79
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# associated salmon output with annotation package
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
names(txi)

#Using DESeq2 to analysis differential gene expression
library(DESeq2)
sampleTable <- data.frame(condition = factor(c(rep("SVZ",2), rep("Control", 8)))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)

# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, name="condition_treated_vs_untreated")
res <- results(dds, contrast=c("condition","treated","untreated"))
res

# Log fold change shrinkage(reduce noises)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC

# p-values and adjusted p-values
resOrdered <- res[order(res$pvalue),] #order our results table by the smallest p value
summary(res)
sum(res$padj < 0.1, na.rm=TRUE) #how many adjusted p-values were less than 0.1
res01 <- results(dds, alpha=0.1) #default alpha is set to 0.1, if the adjusted p value cutoff will be a value other than 0.1, alpha = 0.05
summary(res01)

# MA-plot
plotMA(res, ylim=c(-2,2))

# Plot counts: examine the counts of reads for a single gene across the group 
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# Export results to CVS files
write.csv(as.data.frame(resOrdered), 
          file="desktop/SVZ_control.csv")
          
# Extracting transformed values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

# Effects of transformations on the variance
#(standard deviation of transformed data, across samples, aginst the mean, using the shifted logarithm transformation)library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))
