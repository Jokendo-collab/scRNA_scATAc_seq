library(Signac)
library(Seurat)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)



counts <- Read10X_h5("/Users/jimenezea/OneDrive - National Institutes of Health/Hui_scRNA/24hpa_filtered_feature_bc_matrix.h5")
fragpath <- "/Users/jimenezea/OneDrive - National Institutes of Health/Hui_scRNA/24hpa_atac_fragments.tsv.gz"



pbmc <- CreateSeuratObject(
counts = counts$`Gene Expression`,
assay = "RNA"
)

#Load the zebrafish annotation
gtf <- rtracklayer::import('/Users/jimenezea/OneDrive - National Institutes of Health/Desktop/zebrafish_annotations/chrDanio_rerio.GRCz11.96.gtf')
annotations <- gtf[gtf$type == 'gene']
annotations <- keepStandardChromosomes(annotations, pruning.mode = 'coarse')
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(pbmc) <- annotations OR Annotation <- annotations

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations)
  
  
May have to eventually add chr: 
rownames(peaks) = paste0("chr",rownames(peaks))

  
Error fragments not filtered:
Briefly, run the following cmds in linux/unix terminal :
	1.	gzip -d <fragment.tsv.gz>
	2.	bgzip <fragment.tsv>
	3.	tabix -p bed <fragment.tsv.gz>

To upload zebrafish annotation:
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
gtf <- rtracklayer::import('/Users/jimenezea/OneDrive - National Institutes of Health/Desktop/zebrafish_annotations/chrDanio_rerio.GRCz11.96.gtf')
annotations <- gtf[gtf$type == 'gene']
annotations <- keepStandardChromosomes(annotations, pruning.mode = 'coarse')
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(pbmc) <- annotations



