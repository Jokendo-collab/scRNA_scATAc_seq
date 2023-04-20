library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


pbmc.data <- Read10X(data.dir = "/Users/jimenezea/OneDrive - National Institutes of Health/Hui_scRNA/24hpa_filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 12, cells = 500, balanced = TRUE)
DimHeatmap(pbmc)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap", pt.size = 1)
saveRDS(pbmc, file = "/Users/jimenezea/OneDrive - National Institutes of Health/Hui_scRNA/24h_rna.rds")
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
group_by(cluster) %>%
slice_max(n = 2, order_by = avg_log2FC)
write.table(pbmc.markers, file="/Users/jimenezea/OneDrive - National Institutes of Health/Hui_scRNA/24h_rna_markers.txt", sep = "\t", quote = FALSE)




new.cluster.ids <- c("0", "cl1", "2", "3", "neural tube", "neuromast", "seven", "eight", "nine", "ten", "eleven" )
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

