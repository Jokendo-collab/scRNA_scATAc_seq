#set the working directory
setwd("/Users/okendojo/Desktop/scRNA/pbmcTutorial")
# The script is adapted from: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
#Load the required libraries
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
options(future.globals.maxSize = 1e+09)

#SET UP THE SEURAT OBJECT
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/okendojo/Desktop/scRNA/hui_data /24hpa_outs_040323/24hpa_filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`, project = "notchord", min.cells = 3, min.features = 200)
pbmc

# Lets examine a few genes in the first thirty cells
pbmc.data$`Gene Expression`[c("stm", "col2a1a", "col9a2"), 1:30]

#Standard pre-processing procedure; selecting cells for further analysis
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalizing the data; After removing unwanted cells from the dataset, the next step is to normalize the data. 
pbmc = NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
#We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset 

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# examine and visualize PCA in a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 10)
DimPlot(pbmc)
VizDimLoadings(pbmc)

#Helps in determining the primary source of heterogeneity
DimHeatmap(pbmc,dims = 1:5, cells = 500, balanced = TRUE,nfeatures = 25,fast = TRUE)

#check PC1 - 15
DimHeatmap(pbmc,dims = 1:5, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

#The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line).
JackStrawPlot(pbmc, dims = 1:15)

#elbowplot shows which PC captures the majority of the true signals
ElbowPlot(pbmc)

#cluster the cells

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

#Save the file to avoid reruning the script next time
saveRDS(pbmc, file = "pbmc_tutorial.rds")
list.files()

#Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 20)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#write.csv(pbmc.markers,"notoMarkers.csv")

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#Visualize marker expression
VlnPlot(pbmc, features = c("robo3", "ecrg4b"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations

FeaturePlot(pbmc, features = c("robo3", "ecrg4b"))

#Heatmap plotting
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) 

  # Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) 













































