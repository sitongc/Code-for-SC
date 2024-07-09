a <- read.csv('/Users/chensitong/Desktop/sc/p5-6/p5-6.csv',header = T, row.names = 1)
b <- a[!duplicated(a$GENE, fromLast=TRUE), ] 
write.csv(b,'desktop/p1-3.csv')
bb <- read.csv('desktop/p1-3.csv',row.names = 1,header = T)
cbmc.rna <- as.sparse(a)

##seurat tutorial 
pbmc <- CreateSeuratObject(counts = cbmc.rna, project = 'dpi25', min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#expression level data without mt number 
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

##Normalizing the data
pbmc <- NormalizeData(pbmc)

#Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

##PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

##determine the significant PC
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.61)
pbmc <- RunTSNE(pbmc, dims = 1:11, method = "FIt-SNE")

# Find the markers that define each cluster, and use these to annotate the clusters, we use
# max.cells.per.ident to speed up the process
pbmc.rna.markers <- FindAllMarkers(pbmc, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)

FeaturePlot(pbmc, features = c( 'Nkain4','Naaa','Aqp4','Gja1','Slc1a2','Slc1a3','Olig2','Olig1'), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
DimPlot(pbmc) 