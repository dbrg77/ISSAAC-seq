#====== Author ======#
#=== Weilong Yang ===#

# Load packages
library(Matrix)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
set.seed(1234)
library(monocle)
library(mgsub)
library(stringi)

#===============construct object

# Read RNA data
mat_rep1_RNA <- Read10X(data.dir = "./rep1/RNA/outs/Solo.out/Gene/filtered/")
mat_rep2_RNA <- Read10X(data.dir = "./rep2/RNA/outs/Solo.out/Gene/filtered/")
clonames_rep1_RNA_new <- paste(stri_reverse(mgsub(colnames(mat_rep1_RNA),c("A","T","C","G"), c("T","A","G","C"))), sep = "", "-1")
clonames_rep2_RNA_new <- paste(stri_reverse(mgsub(colnames(mat_rep2_RNA),c("A","T","C","G"), c("T","A","G","C"))), sep = "", "-2")
colnames(mat_rep1_RNA) <- clonames_rep1_RNA_new
colnames(mat_rep2_RNA) <- clonames_rep2_RNA_new
mat_all_RNA <- cbind(mat_rep1_RNA, mat_rep2_RNA[rownames(mat_rep1_RNA),])

# Read ATAC data
mat_all_ATAC <- Read10X_h5("./mCortex_ATAC_aggr/outs/filtered_peak_bc_matrix.h5")

# Find barcodes both in RNA and ATAC data
RNA_BC <- colnames(mat_all_RNA)
ATAC_BC <- colnames(mat_all_ATAC)
cor_BC <- intersect(RNA_BC,ATAC_BC)
mat_cor_all_RNA <- mat_all_RNA[,cor_BC]
mat_cor_all_ATAC <- mat_all_ATAC[,cor_BC]

# Create Seurat object
mCortex_all <- CreateSeuratObject(mat_cor_all_RNA, project = "mCortex")
mCortex_all[["percent.mt"]] <- PercentageFeatureSet(mCortex_all, pattern = "^mt-")

# Add ATAC-seq data
grange.counts <- StringToGRanges(rownames(mat_cor_all_ATAC), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
mat_cor_all_ATAC <- mat_cor_all_ATAC[as.vector(grange.use),]

# Add gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotations) <- ucsc.levels
genome(annotations) <- "mm10"
frag.file <- "./mCortex_ATAC_aggr/outs/fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = mat_cor_all_ATAC,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  annotation = annotations
)
mCortex_all[["ATAC"]] <- chrom_assay
DefaultAssay(mCortex_all) <- "ATAC"
mCortex_all <- TSSEnrichment(mCortex_all)
mCortex_all <- NucleosomeSignal(mCortex_all)
mCortex_all$blacklist_fraction <- FractionCountsInRegion(
  object = mCortex_all,
  assay = 'ATAC',
  regions = blacklist_mm10
)

# Add group information
Group <- c(rep("mCortex_rep1",3974),rep("mCortex_rep2",6404))
mCortex_all <- AddMetaData(mCortex_all, Group, col.name = "Group")
saveRDS(mCortex_all, "./analysis_pipe_out/mCortex_all_original.rds")

# Export expression matrix
expression_matrix <- mCortex_all@assays$RNA@counts
write.csv(expression_matrix, "./analysis_pipe_out/mCortex_all_expression_matrix.csv", row.names = T, col.names = T, quote=F)

#===============RNA analysis

mCortex_all <- readRDS("./analysis_pipe_out/mCortex_all_original.rds")
DefaultAssay(mCortex_all) <- "RNA"
mCortex_all <- FindVariableFeatures(mCortex_all, nfeatures = 3000)
mCortex_all <- NormalizeData(mCortex_all)
mCortex_all <- ScaleData(mCortex_all)
mCortex_all <- RunPCA(mCortex_all, npcs = 30)
mCortex_all <- RunUMAP(mCortex_all, dims = 1:30, reduction.name = "umap.rna")
mCortex_all <- FindNeighbors(mCortex_all, dims = 1:30, k.param = 30)
mCortex_all <- FindClusters(mCortex_all, resolution = 1.0, algorithm = 3)
DimPlot(mCortex_all, reduction = "umap.rna", label = TRUE, group.by = "RNA_snn_res.1")  + ggtitle("RNA UMAP")

# Export the RNA UMAP coordinates
write.csv(mCortex_all@reductions$umap.rna@cell.embeddings, "./analysis_pipe_out/mCortex_all_RNA_UMAP_coordinates.csv", row.names = T, col.names = T, quote=F)

# Find different genes
Idents(mCortex_all) <- "RNA_snn_res.1"
diff.wilcox = FindAllMarkers(mCortex_all)
write.csv(diff.wilcox, "./analysis_pipe_out/mCortex_all_diff_genes_wilcox_r1.0.csv", row.names = F)

#===============ATAC analysis

DefaultAssay(mCortex_all) <- "ATAC"
mCortex_all <- RunTFIDF(mCortex_all)
mCortex_all <- FindTopFeatures(mCortex_all, min.cutoff = 'q0')
mCortex_all <- RunSVD(object = mCortex_all)
mCortex_all <- RunUMAP(
  object = mCortex_all,
  reduction = 'lsi',
  dims = 2:30, reduction.name = "umap.atac"
)
mCortex_all <- FindNeighbors(
  object = mCortex_all,
  reduction = 'lsi',
  dims = 2:30, k.param = 30
)
mCortex_all <- FindClusters(
  object = mCortex_all,
  resolution = 1.0, algorithm = 3
)
DimPlot(mCortex_all, reduction = "umap.atac", label = TRUE, group.by = "ATAC_snn_res.1") + ggtitle("ATAC UMAP")

#Export the ATAC UMAP coordinates
write.csv(mCortex_all@reductions$umap.atac@cell.embeddings, "./analysis_pipe_out/mCortex_all_ATAC_UMAP_coordinates.csv", row.names = T, col.names = T, quote=F)

#Find different peaks
diff.LR_peaks = FindAllMarkers(mCortex_all, assay = 'ATAC',test.use = "LR")
write.csv(diff.LR_peaks, "./analysis_pipe_out/mCortex_all_diff_peaks_LR_r1.0.csv", row.names = F)

#===============Annotate cell types in the dataset by transferring labels from an existing mCortex-seq dataset for the adult mouse mCortex, produced by the Allen Institute.
## Note the "allen_mCortex.rds" file comes from the R package Signac (https://github.com/timoast/signac)
## See here: https://github.com/timoast/signac/blob/master/vignettes/mouse_brain_vignette.Rmd

allen <- readRDS("./analysis_pipe_out/allen_mCortex.rds")

#use the RNA assay in the mCortex-seq data for integration with mCortex-seq
DefaultAssay(mCortex_all) <- 'RNA'
transfer.anchors <- FindTransferAnchors(
  reference = allen,
  query = mCortex_all,
  dims = 1:30,
  reduction = 'cca'
)
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen$subclass,
  weight.reduction = mCortex_all[['pca']],
  dims = 1:30
)
mCortex_all <- AddMetaData(object = mCortex_all, metadata = predicted.labels)
saveRDS(mCortex_all, "./analysis_pipe_out/mCortex_all.final.rds")

# Export metadata
write.csv(mCortex_all@meta.data, "./analysis_pipe_out/mCortex_all_metadata.csv", quote = FALSE)

#===============Pseudotime analysis (RNA cluster R14 and R20)

#RNA analysis

## Pick up RNA cluster RNA_snn_res.1 == 14 and RNA_snn_res.1 == 20
mCortex_all <- readRDS("./analysis_pipe_out/mCortex_all.final.rds")
DefaultAssay(mCortex_all) <- "RNA"
mCortex_all_R14_R20_id <- subset(mCortex_all@meta.data, RNA_snn_res.1 %in% c(14, 20))
mCortex_all_R14_R20 <- subset(mCortex_all, cells=row.names(mCortex_all_R14_R20_id))
DimPlot(mCortex_all_R14_R20, group.by = "RNA_snn_res.1", label = TRUE, reduction = "umap.rna") + NoLegend() + ggtitle("RNA UMAP")
DefaultAssay(mCortex_all_R14_R20) <- "RNA"
Idents(mCortex_all_R14_R20) <- "RNA_snn_res.1"

# Construct monocle object
data <- as(as.matrix(mCortex_all_R14_R20@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = mCortex_all_R14_R20@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds_RNA <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
mycds_RNA <- estimateSizeFactors(mycds_RNA)
mycds_RNA <- estimateDispersions(mycds_RNA, cores=4, relative_expr = TRUE)
diff.wilcox_mCortex_all_R14_R20 = FindMarkers(mCortex_all_R14_R20, ident.1 = "14", ident.2 = "20") # Find different genes between RNA clusters R14 and R20
write.csv(diff.wilcox_mCortex_all_R14_R20, "./analysis_pipe_out/diff.wilcox_mCortex_all_R14_R20.csv")
diff.genes <- subset(diff.wilcox_mCortex_all_R14_R20, p_val_adj < 0.05) %>% rownames()
mycds_RNA <- setOrderingFilter(mycds, diff.genes)
mycds_RNA <- reduceDimension(mycds_RNA, method = 'DDRTree', max_components = 2)
mycds_RNA <- orderCells(mycds_RNA)
p1 <- plot_cell_trajectory(mycds_RNA, color_by = "State", cell_size = 0.1)
p2 <- plot_cell_trajectory(mycds_RNA, color_by = "RNA_snn_res.1", cell_size = 0.1)
p3 <- plot_cell_trajectory(mycds_RNA, color_by = "Pseudotime", cell_size = 0.1)
p1 | p2 | p3

## Export the RNA peudotime information
peudotime_inf = as.data.frame(t(mycds@reducedDimS))
colnames(peudotime_inf) = c("component1", "component2")
peudotime_inf$Pseudotime = mycds$Pseudotime
peudotime_inf$State = mycds$State
write.csv(peudotime_inf, "./analysis_pipe_out/mCortex_all_OPC_to_Oligo_peudotime_information.csv")

# ATAC analysis

## Create a gene activity score matrix
DefaultAssay(mCortex_all_R14_R20) <- "ATAC"
gene.activities <- GeneActivity(mCortex_all_R14_R20)
mCortex_all_R14_R20[['activity']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(mCortex_all_R14_R20) <- "activity"
mCortex_all_R14_R20 <- NormalizeData(mCortex_all_R14_R20)
data <- as(as.matrix(mCortex_all_R14_R20@assays$activity@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = mCortex_all_R14_R20@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds_ATAC <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
mycds_ATAC <- estimateSizeFactors(mycds_ATAC)
mycds_ATAC <- estimateDispersions(mycds_ATAC, cores=4, relative_expr = TRUE)
mycds_ATAC <- setOrderingFilter(mycds_ATAC, diff.genes)
mycds_ATAC <- reduceDimension(mycds_ATAC, max_components = 2, method = 'DDRTree')
mycds_ATAC <- orderCells(mycds_ATAC, reverse=T)
p4 <- plot_cell_trajectory(mycds_ATAC, color_by = "State", cell_size = 0.1)
p5 <- plot_cell_trajectory(mycds_ATAC, color_by = "RNA_snn_res.1", cell_size = 0.1)
p6 <- plot_cell_trajectory(mycds_ATAC, color_by = "Pseudotime", cell_size = 0.1)
p4 | p5 | p6

# Find genes both in RNA expression matrix and gene activity score matrix
RNA_gene_list <- rownames(mycds_RNA)
activity_gene_list <- rownames(mycds_ATAC)
cor_gene <- intersect(RNA_gene_list,activity_gene_list)
diff.genes_2 <- subset(diff.wilcox_mCortex_all_R14_R20,p_val_adj < 0.1) %>% rownames()
final.gene <- intersect(cor_gene,diff.genes_2)

# Analysis gene expression and gene activity score change with RNA pseudotime

## Analysis gene expression change with RNA pseudotime
newdata <- data.frame(Pseudotime = seq(min(pData(mycds_RNA)$Pseudotime),
                                       max(pData(mycds_RNA)$Pseudotime), length.out = 100)) #normal RNA pseudotime
genSmoothCurves_mat_RNA <- genSmoothCurves(mycds[final.gene,],
                                           new_data = newdata,
                                           cores = 10) #use RNA normal pseudotime to smooth gene expression
write.csv(genSmoothCurves_mat_RNA, "./analysis_pipe_out/mCortex_all_OPC_to_Oligo_cor_gene_expression_pseudotime.csv")

## Analysis gene activity score change with RNA pseudotime
genSmoothCurves_mat_activity <- genSmoothCurves(mycds_ATAC[final.gene,],
                                                new_data = newdata,
                                                cores = 10) # use RNA normal pseudotime to smooth gene activity score
write.csv(genSmoothCurves_mat_activity, "./analysis_pipe_out/mCortex_all_OPC_to_Oligo_cor_gene_activity_score_pseudotime.csv")
