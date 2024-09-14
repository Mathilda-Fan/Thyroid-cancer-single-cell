#Pseudotime trajectory analysis
library(Seurat)
library(monocle)
library(igraph)

load("scRNA_harmony.Rdata")
data <- as(as.matrix(scRNA_harmony_V3@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA_harmony_V3@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))
HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)
plot1 <- plot_cell_trajectory(HSMM, color_by = "Annotation1")
plot1 <- plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot1 <- plot_cell_trajectory(HSMM, color_by = "grade")
plot1 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")

HSMM6 <- ds
plot1 <- plot_cell_trajectory(HSMM, color_by = "group") +
  facet_wrap(~group, nrow = 1)
plot1
ggsave("trajectory_group_row.pdf", plot = plot1,width = 10,height = 7)
plot1 <- plot_cell_trajectory(HSMM, color_by = "State") +
  facet_wrap(~State, nrow = 1)

blast_genes <- row.names(subset(fData(HSMM),
                                gene_short_name %in% c("ANXA2","SOX4")))
plot1 <- plot_genes_jitter(HSMM[blast_genes,],
                  grouping = "grade",
                  min_expr = 0.1)

HSMM_expressed_genes <-  row.names(subset(fData(HSMM),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("ANXA2")))
cds_subset <- HSMM_filtered[my_genes,]
plot1 <- plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")
plot1 <- plot_genes_in_pseudotime(cds_subset, color_by = "Annotation1")
plot1 <- plot_genes_in_pseudotime(cds_subset, color_by =  "grade")