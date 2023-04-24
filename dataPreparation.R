
setwd("/Users/okendojo/Desktop/scRNA/hui_data /24hpa_outs_040323")
library(Seurat)
library(reticulate)

# load the data
counts <- Read10X_h5("24hpa_filtered_feature_bc_matrix.h5")

rownames(counts$`Gene Expression`) <- make.unique(rownames(counts$`Gene Expression`))

rna <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = 'RNA', min.cells = 5, min.features = 500, project = 'scRNAseq')

scrublet <- read.table("24hpa_filtered_feature_bc_matrix/barcodes.tsv.gz", sep = "\t")
rownames(scrublet) <- colnames(counts)
rna <- AddMetaData(rna, metadata = scrublet)
rna <- RenameCells(rna, add.cell.id = 'rna')
mito.features <- grep(pattern = "^MT-", x = rownames(x = rna), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = rna, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = rna, slot = 'counts'))
rna$percent.mito <- percent.mito

# run the quality analysis
rna <- subset(x = rna, subset = nCount_RNA > 200 & nCount_RNA < 20000 & percent.mito < 0.2)


# preprocessing
rna = NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

#rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures = 3000)
rna <- ScaleData(rna)
rna <- RunPCA(rna, npcs = 50)
rna <- RunTSNE(rna, dims = 1:30)
rna <- FindNeighbors(rna, dims = 1:30)
rna <- FindClusters(rna, resolution = 0.4, algorithm = 3)
rna <- RunUMAP(rna, dims = 1:10)

#Add cell types to the markers
new.cluster.ids <- c(
  "Blastema cells",
  'Epidermal cells',
  'Granulocytes',
  'Proliferating cells',
  'Double negative T cell', ####
  'Keratenocyte',
  'Blastema',
  'CD8 effector', ####
  'Osteoblast',
  'CD16+ Monocytes', ####
  'Dendritic cell', ####
  'pDC',####
  'Platelet')####

names(x = new.cluster.ids) <- levels(x = rna)
rna <- RenameIdents(object = rna, new.cluster.ids)
rna$celltype <- Idents(rna)
#nk.cells <- subset(rna, subset = celltype == 'CD8 effector')
#gzmk <- GetAssayData(nk.cells, assay = 'RNA', slot = 'data')['stm', ]
#nk.cells$bright <- ifelse(gzmk > 1, 'NK bright', 'NK dim')
ctypes <- as.vector(rna$celltype)
names(ctypes) <- names(rna$celltype)
#ctypes[Cells(nk.cells)] <- nk.cells$bright
rna <- AddMetaData(rna, metadata = ctypes, col.name = 'celltype')
saveRDS(rna, "pbmc_10k_v3.rds")







