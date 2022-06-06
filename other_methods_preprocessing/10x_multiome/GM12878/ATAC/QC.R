#! /usr/bin/env Rscript

# Load packages
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)

# Input data
input_rep1 <- Read10X_h5(data.dir = "GSM5085810_GM12878_rep1_filtered_feature_bc_matrix.h5")
input_rep2 <- Read10X_h5(data.dir = "GSM5085812_GM12878_rep2_filtered_feature_bc_matrix.h5")
rna_counts_rep1 <- input_rep1$`Gene Expression`
atac_counts_rep1 <- input_rep1$Peaks
rna_counts_rep2 <- input_rep2$`Gene Expression`
atac_counts_rep2 <- input_rep2$Peaks

# Create Seurat object
GM12878_rep1 <- CreateSeuratObject(counts = rna_counts_rep1)
GM12878_rep1[["percent.mt"]] <- PercentageFeatureSet(GM12878_rep1, pattern = "^MT-")
GM12878_rep2 <- CreateSeuratObject(counts = rna_counts_rep2)
GM12878_rep2[["percent.mt"]] <- PercentageFeatureSet(GM12878_rep2, pattern = "^MT-")

# Add ATAC-seq data
grange.conuts_rep1 <- StringToGRanges(rownames(atac_counts_rep1), sep = c(":", "-"))
grange.use_rep1 <- seqnames(grange.conuts_rep1) %in% standardChromosomes(grange.conuts_rep1)
atac_counts_rep1 <- atac_counts_rep1[as.vector(grange.use_rep1),]
grange.conuts_rep2 <- StringToGRanges(rownames(atac_counts_rep2), sep = c(":", "-"))
grange.use_rep2 <- seqnames(grange.conuts_rep2) %in% standardChromosomes(grange.conuts_rep2)
atac_counts_rep2 <- atac_counts_rep2[as.vector(grange.use_rep2),]

# Add gene annotations 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotations) <- ucsc.levels
genome(annotations) <- "hg38"
frag.file_rep1 <- "/path/GSM5085811_GM12878_rep1_atac_fragments.tsv.gz"
frag.file_rep2 <- "/path/GSM5085811_GM12878_rep2_atac_fragments.tsv.gz"
chrom_assay_rep1 <- CreateChromatinAssay(
  counts = atac_counts_rep1,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file_rep1,
  annotation = annotations
)
GM12878_rep1[["ATAC"]] <- chrom_assay_rep1
write.csv(GM12878_rep1@meta.data, "GM12878_rep1_metadata.csv", row.names = T, col.names = T, quote=F)
chrom_assay_rep2 <- CreateChromatinAssay(
  counts = atac_counts_rep2,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file_rep2,
  annotation = annotations
)
GM12878_rep2[["ATAC"]] <- chrom_assay_rep2
write.csv(GM12878_rep2@meta.data, "GM12878_rep2_metadata.csv", row.names = T, col.names = T, quote=F)




