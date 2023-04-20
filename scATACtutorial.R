#set working directory
#this workflow was adapted from: https://stuartlab.org/signac/articles/pbmc_vignette.html
setwd("/Users/okendojo/Desktop/scRNA/hui_data /24hpa_outs_040323")

#Load all the libraries needed for the analysis
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(org.Dr.eg.db)
set.seed(1234)


#PRE-PROCESSING WORKFLOW
counts <- Read10X_h5(filename = "24hpa_filtered_feature_bc_matrix.h5")
metadata <- read.csv(
  file = "24hpa_per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)


chrom_assay <- CreateChromatinAssay( counts = counts$Peaks, sep = c(":", "-"),
                                     #genome = 'danRer11', 
                                     fragments = '24hpa_atac_fragments.tsv.gz',
                                     min.cells = 10, min.features = 200, )


pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

#Check the object
pbmc
pbmc[['peaks']]

#see the genomic ranges associated with each feature in the object
granges(pbmc)

#We can also add gene annotations to the pbmc object for the human genome. 
#This will allow downstream functions to pull the gene annotation information directly from the object.
# extract gene annotations from EnsDb
gtf <- rtracklayer::import('Danio_rerio.GRCz11.109.gtf')
annotations <- gtf[gtf$type == 'gene']
annotations <- keepStandardChromosomes(annotations, pruning.mode = 'coarse')
seqlevelsStyle(annotations) <- 'NCBI'
Annotation(pbmc) <- annotations

head(Annotation(pbmc))
head(Fragments(pbmc)[[1]])


#Computing QC metrices
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)


# # add blacklist ratio and fraction of reads in peaks
# pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
# pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

#Inspect the enrichment score by plotting
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

#We can also look at the fragment length periodicity for all the cells, and group by cells with high or low nucleosomal signal strength.
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

#Finally we remove cells that are outliers for these QC metrics.
pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
pbmc

#Normalization and linear dimensional reduction
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

#Visualize the SVD
DepthCor(pbmc)

#Non-linear dimension reduction and clustering
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) 

#Create a gene activity matrix
gene.activities <- GeneActivity(pbmc)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

#Now we can visualize the activities of canonical marker genes to help interpret our ATAC-seq clusters
DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

#Integrating with scRNA-seq data
# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("pbmc_10k_v3.rds")

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

#Visualize the integration
plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

#We note that cluster 14 maps to CD4 Memory T cells, but is a very small cluster with lower QC metrics. 
#As this group is likely representing low-quality cells, we remove it from downstream analysis.

pbmc <- subset(pbmc, idents = 14, invert = TRUE)
pbmc <- RenameIdents(
  object = pbmc,
  '0' = 'CD14 Mono',
  '1' = 'CD4 Memory',
  '2' = 'CD8 Effector',
  '3' = 'CD4 Naive',
  '4' = 'CD14 Mono',
  '5' = 'DN T',
  '6' = 'CD8 Naive',
  '7' = 'NK CD56Dim',
  '8' = 'pre-B',
  '9' = 'CD16 Mono',
  '10' = 'pro-B',
  '11' = 'DC',
  '12' = 'NK CD56bright',
  '13' = 'pDC'
)

#Find differentially accessible peaks between clusters
# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14 Mono",
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)

#visualize differential accessible peaks
plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14 Mono")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

#Another way to find DA regions between two groups of cells is to look at the fold change accessibility between two groups of cells.

fc <- FoldChange(pbmc, ident.1 = "CD4 Naive", ident.2 = "CD14 Mono")
head(fc)

#Plotting genomic regions
#We can plot the frequency of Tn5 integration across regions of the genome for cells grouped by cluster, cell type, or
#any other metadata stored in the object for any genomic region using the CoveragePlot() function. 

# set plotting order
levels(pbmc) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 Effector","DN T","NK CD56bright","NK CD56Dim","pre-B",'pro-B',"pDC","DC","CD14 Mono",'CD16 Mono')

#genonome coverage
CoveragePlot(
  object = pbmc,
  region = rownames(da_peaks)[1:3],
  extend.upstream = 40000,
  extend.downstream = 20000
)


BiocManager::install("EnsDb.danRer11")








































































































































