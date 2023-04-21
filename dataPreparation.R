
setwd("/Users/okendojo/Desktop/scRNA/hui_data /24hpa_outs_040323")
library(Seurat)
library(reticulate)

# load 10x v3 pbmc data
counts <- Read10X_h5("24hpa_filtered_feature_bc_matrix.h5")

rownames(counts$`Gene Expression`) <- make.unique(rownames(counts$`Gene Expression`))

rna <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = 'RNA', min.cells = 5, min.features = 500, project = '10x_RNA')

scrublet <- read.table("24hpa_filtered_feature_bc_matrix/barcodes.tsv.gz", sep = "\t")
rownames(scrublet) <- colnames(counts)
rna <- AddMetaData(rna, metadata = scrublet)
rna <- RenameCells(rna, add.cell.id = 'rna')
mito.features <- grep(pattern = "^MT-", x = rownames(x = rna), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = rna, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = rna, slot = 'counts'))
rna$percent.mito <- percent.mito

# QC
rna <- subset(x = rna, subset = nCount_RNA > 2000 & nCount_RNA < 20000 & percent.mito < 0.2)# & observed < 0.1)

# preprocessing
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures = 3000)
rna <- ScaleData(rna)
rna <- RunPCA(rna, npcs = 100)
rna <- RunTSNE(rna, dims = 1:30)
rna <- FindNeighbors(rna, dims = 1:30)
rna <- FindClusters(rna, resolution = 0.4, algorithm = 3)
rna <- RunUMAP(rna, dims = 1:10)

#rna <- RunUMAP(object = rna, graph = 'RNA_nn', metric = 'euclidean')



new.cluster.ids <- c(
  "CD14+ Monocytes",
  'CD4 Memory',
  'CD4 Naive',
  'pre-B cell',
  'Double negative T cell',
  'NK cell',
  'B cell progenitor',
  'CD8 effector',
  'CD8 Naive',
  'CD16+ Monocytes',
  'Dendritic cell', 
  'pDC',
  'Platelet'
)

names(x = new.cluster.ids) <- levels(x = rna)
rna <- RenameIdents(object = rna, new.cluster.ids)
rna$celltype <- Idents(rna)
nk.cells <- subset(rna, subset = celltype == 'NK cell')
gzmk <- GetAssayData(nk.cells, assay = 'RNA', slot = 'data')['stm', ]
nk.cells$bright <- ifelse(gzmk > 1, 'NK bright', 'NK dim')
ctypes <- as.vector(rna$celltype)
names(ctypes) <- names(rna$celltype)
ctypes[Cells(nk.cells)] <- nk.cells$bright
rna <- AddMetaData(rna, metadata = ctypes, col.name = 'celltype')
saveRDS(rna, "pbmc_10k_v3.rds")





