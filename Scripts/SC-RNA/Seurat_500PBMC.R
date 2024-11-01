# load seurat package
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(dplyr)
library(SingleCellExperiment)


# load cell ranger output 
data_dir <- 'desktop/filtered_feature_bc_matrix/'
pbmc <- Read10X(data.dir = data_dir)

# create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc, project = "500_pbmc")

# QC based on the MT
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# keep the RNA feature between 200 to 7500 and percentige of MT less than 10
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Normaling the data by log 2 folder change
pbmc <- NormalizeData(pbmc)

# Identification of highly variable features nfeatures > 7000.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5500)

top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

#determine how many pc do we use
ElbowPlot(pbmc)

# cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

# load the reference from celldex
reference <- celldex::HumanPrimaryCellAtlasData()

# Convert seurat object to single cell experiment object
sce <- as.SingleCellExperiment(pbmc)

singleR_labels <- SingleR(test = sce, ref = reference, labels = reference$label.main)
pbmc$SingleR_label <- singleR_labels$labels

DimPlot(pbmc, reduction = "umap", group.by = "SingleR_label", label = TRUE, label.size = 4) +
  ggtitle("500 pbmc")

# create the heatmap based on the SingleR annotation.
Idents(pbmc) <- pbmc$SingleR_label

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
