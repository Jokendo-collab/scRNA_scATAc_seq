#Set the working directory
setwd("/Users/okendojo/Desktop/scRNA/hui_data /reprocessedData")

#Load the libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(destiny)
library(presto)
library(loomR)
library(reticulate)

# 1. Create the seurat object
counts <- Read10X(data.dir = "36hpa_filtered_feature_bc_matrix/")
seurat <- CreateSeuratObject(counts$`Gene Expression`, project="36hpa")

# 2.Quality control; It is done to remove MT DNA, cells with too few and too many detected genes
#Calculate mt transcript %
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")

#Look at the distribution of the features
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # Shows the dotplots
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0) # no dotplots

#Detect the number of detected genes and features
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seurat <- subset(seurat, subset = nFeature_RNA > 5 & nFeature_RNA < 4000 & percent.mt < 2)

# 3. Normalize the data
seurat <- NormalizeData(seurat)

# 4. Step 4. Feature selection for following heterogeneity analysis
seurat <- FindVariableFeatures(seurat, nfeatures = 2000)

#visualize the result in a variable feature plot
top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2

# 5. Data scaling
seurat <- ScaleData(seurat)

# Step 6. Linear dimensionality reduction using principal component analysis (PCA)
seurat <- RunPCA(seurat, npcs = 50)
DimPlot(seurat)
DimHeatmap(seurat)

#Visualize the informative PCs
ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))

#check the distribution of genes in the top PCs
PCHeatmap(seurat, dims = 1:5, cells = 500, balanced = TRUE, ncol = 4)

#Step 7. Non-linear dimension reduction for visualization
seurat <- RunTSNE(seurat, dims = 1:20)
seurat <- RunUMAP(seurat, dims = 1:20)

#Visualize the results
plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)
plot1 + plot2

#check whether certain cell types or cell states exist in the data
plot1 <- FeaturePlot(seurat, c("col2a1a","stm","col9a2","mgp","krt94","adgrg6",
                               "ecrg4b",
                               "itga8"),ncol=3, reduction = "tsne")
plot2 <- FeaturePlot(seurat, c("col2a1a","stm","col9a2","mgp","krt94","adgrg6",
                               "ecrg4b",
                               "itga8"),ncol=3, reduction = "umap")

plot1 / plot2

# Step 8. Cluster the cells
seurat <- FindNeighbors(seurat, dims = 1:20)

seurat <- FindClusters(seurat, resolution = 1)

#Visualize the clusters
plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE) + NoLegend()
plot2 <- DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()
plot1 + plot2

# Step 9. Annotate cell clusters

ct_markers = c("col2a1a","stm","col9a2","mgp","krt94","adgrg6",
               "ecrg4b","pdlim3b","pou2f2a.1","col14a1a","mgst1.2","evpla","FP102018.1","ctsk","trpv1","plpp3",
               "itga8")
DoHeatmap(seurat, features = ct_markers) + NoLegend()

#identify cluster markers for each of the cell cluster identified.
cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
cl_markers %>% group_by(cluster) %>% top_n(n = 2)

cl_markers_presto <- wilcoxauc(seurat)
cl_markers_presto %>%
  filter(logFC > log(1.5) & pct_in > 20 & padj < 0.05) %>%
  group_by(group) %>%
  arrange(desc(logFC), .by_group=T) %>%
  top_n(n = 2, wt = logFC) %>%
  print(n = 40, width = Inf)
#write.csv(cl_markers_presto,"24hpa.csv")

#No matter with which method, the identified top cluster markers can be next visualized by a heatmap

top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10)#, wt = avg_logFC)
DoHeatmap(seurat, features = top10_cl_markers$gene) + NoLegend()

#check those markers of different clusters in more details, by doing feature plot or violin plot 

plot1 <- FeaturePlot(seurat, c("robo3","col5a1"), ncol = 1)
plot2 <- VlnPlot(seurat, features = c("robo3","col5a1"), pt.size = 0)
plot1 + plot2 + plot_layout(widths = c(1, 2))

#We can replace the cell cluster labels by the annotation, but this is optional
new_ident <- c("Dorsal telen. NPC",
                        "Midbrain-hindbrain boundary neuron",
                        "Dorsal telen. neuron",
                        "Dien. and midbrain excitatory neuron",
                        "MGE-like neuron","G2M dorsal telen. NPC",
                        "Dorsal telen. IP","Dien. and midbrain NPC",
                        "Dien. and midbrain IP and excitatory early neuron",
                        "G2M Dien. and midbrain NPC",
                        "G2M dorsal telen. NPC",
                        "Dien. and midbrain inhibitory neuron",
                        "Dien. and midbrain IP and early inhibitory neuron",
                        "Ventral telen. neuron",
                        "Unknown 1",
                        "Unknown 2",
                        "Unknown 3",
                        "Unknown 4",
                        "Unknown 5",
                        "Dorsal telen. NPC",
                        "Midbrain-hindbrain boundary neuron")

names(x = new_ident) <- levels(x = seurat)

seurat <- RenameIdents(object = seurat, new_ident)
seurat$celltype <- Idents(seurat)
ctypes <- as.vector(seurat$celltype)
names(ctypes) <- names(seurat$celltype)
seurat$celltype <- AddMetaData(seurat, metadata = ctypes, col.name = 'celltype')
saveRDS(seurat, "24hpa.rds")


DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()


#==========================================================

#Step 10. Pseudotemporal cell ordering
#First of all, cells of interest are extracted. Afterwards, we re-identify highly variable genes for the subset cells,
#as genes representing differences between dorsal telencephalic cells and other cells are no longer informative
seurat_dorsal <- subset(seurat, subset = RNA_snn_res.1 %in% c(0,2,5,6,10))
seurat_dorsal <- FindVariableFeatures(seurat_dorsal, nfeatures = 2000)

#exclude cell-cycle ralated genes since we are interested in molecular mechanism of regeneration
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), unlist(cc.genes))

#We can then check how the data look like, by creating a new UMAP embedding and do some feature plots
seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:20)
FeaturePlot(seurat_dorsal, c("stm","col9a2","mgp","adgtg6"), ncol = 4)

#further reduce the cell cycle influence
seurat_dorsal <- CellCycleScoring(seurat_dorsal,
                                  s.features = cc.genes$s.genes,
                                  g2m.features = cc.genes$g2m.genes,
                                  set.ident = TRUE)
seurat_dorsal <- ScaleData(seurat_dorsal, vars.to.regress = c("S.Score", "G2M.Score"))

seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:20)
FeaturePlot(seurat_dorsal, c("col1a1a","robo3","col5a2a"), ncol = 4)

#Now let's try to run diffusion map to get the cells ordered.

dm <- DiffusionMap(Embeddings(seurat_dorsal, "pca")[,1:20])
dpt <- DPT(dm)
seurat_dorsal$dpt <- rank(dpt$dpt)
FeaturePlot(seurat_dorsal, c("col1a1a","robo3","col5a2a"), ncol=4)

#To visualize expression changes along the constructed pseudotime, a scatter plot with fitted curve is usually a straightforward way.
library(ggplot2)
plot1 <- qplot(seurat_dorsal$dpt, as.numeric(seurat_dorsal@assays$RNA@data["robo3",]),
               xlab="Dpt", ylab="Expression", main="robo3") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot2 <- qplot(seurat_dorsal$dpt, as.numeric(seurat_dorsal@assays$RNA@data["col1a1a",]),
               xlab="Dpt", ylab="Expression", main="col1a1a") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot3 <- qplot(seurat_dorsal$dpt, as.numeric(seurat_dorsal@assays$RNA@data["col5a2a",]),
               xlab="Dpt", ylab="Expression", main="col5a2a") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot1 + plot2 + plot3

#Step 11. Save the result
#saveRDS(seurat, file="UI.rds")
#saveRDS(seurat_dorsal, file="UI_filtered.rds")

# seurat <- readRDS("seurat_obj_all.rds")
# seurat_dorsal <- readRDS("seurat_obj_dorsal.rds")


# Now starts Part 2: when you need to jointly analyze multiple scRNA-seq data sets

#Step 0. Load data
UI = readRDS("UI.rds")
hpa_24 = readRDS("24hpa.rds")
hpa_36 = readRDS("36hpa.rds")
hpa_48 = readRDS("48hpa.rds")

#Step 1. Merge the data sets

seurat <- merge(hpa_48 , y = c(UI,hpa_36,hpa_24)) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)

plot1 <- DimPlot(seurat, group.by="orig.ident")
plot2 <- FeaturePlot(seurat, c("col1a1a","robo3","col5a2a"), ncol=2, pt.size = 0.1)
plot1 + plot2 + plot_layout(widths = c(1.5, 2))

#Step 2-1. Data integration using Seurat
# You should normalize the data if it has not been done the we identify anchors of data sets

seurat_objs <- list(UI = UI,hpa_24 = hpa_24, hpa_36 = hpa_36, hpa_48 = hpa_48)
anchors <- FindIntegrationAnchors(object.list = seurat_objs, dims = 1:30)

#Next, the identified anchor set is passed to the the IntegrateData function to do the expression level correction.
seurat <- IntegrateData(anchors, dims = 1:30)

#Next, we just take the corrected Seurat object and re-run the procedure in Part 1, except for the first two steps (normalization and highly variable gene identification) which should be skipped here

seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 50)
seurat <- RunUMAP(seurat, dims = 1:20)
seurat <- FindNeighbors(seurat, dims = 1:20) %>% FindClusters(resolution = 0.6)

# You may also want to save the object
#saveRDS(seurat, file="integrated_seurat.rds")

# run tSNE/UMAP embedding and clustering 
DefaultAssay(seurat) <- "RNA"
plot1 <- UMAPPlot(seurat, group.by="orig.ident")
plot2 <- UMAPPlot(seurat, label = T)
plot3 <- FeaturePlot(seurat, c("col1a1a","robo3","s100t"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))




































































