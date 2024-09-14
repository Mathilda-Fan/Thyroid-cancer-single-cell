#Spatial transcriptomics sequencing
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

Thyrocytes <- readRDS("D:\\single-wzynew\\Thyrocytes.RDS")
######################
SpatialData <- Load10X_Spatial(data.dir = 'D:\\single-wzynew\\spatial_data\\A1032775DSQ\\',filename = 'filtered_feature_bc_matrix.h5',assay = 'Spatial',slice = 'slice1',filter.matrix = T,to.upper = F,image = NULL)
SpatialData <- SCTransform(SpatialData, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
DefaultAssay(SpatialData) <- "SCT"
saveRDS(SpatialData,'ST-DSQ.rds')
table(Thyrocytes$MERTKEXP)
anchors <- FindTransferAnchors(reference = Thyrocytes, query = SpatialData, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = Thyrocytes$MERTKEXP, prediction.assay = TRUE,
                                  weight.reduction = SpatialData[["pca"]], dims = 1:20)
SpatialData[["predictions"]] <- predictions.assay
DefaultAssay(SpatialData) <- "predictions"
SpatialFeaturePlot(SpatialData, features = c('MERTKpos','MERTKneg'), pt.size.factor =3, ncol = 2, crop = TRUE, alpha = c(0.3, 1))


predictions.assay1 <- TransferData(anchorset = anchors, refdata = Thyrocytes$seurat_clusters, prediction.assay = TRUE,
                                  weight.reduction = SpatialData[["pca"]], dims = 1:20)
SpatialData[["predictions1"]] <- predictions.assay1
DefaultAssay(SpatialData) <- "predictions1"
SpatialFeaturePlot(SpatialData, features = c('0','1','2','3','4','5','6','7'), pt.size.factor =3, ncol = 4, crop = TRUE, alpha = c(0.3, 1))

SpatialData <- RunPCA(SpatialData, assay = "SCT", verbose = FALSE)
SpatialData <- FindNeighbors(SpatialData, reduction = "pca", dims = 1:30)
SpatialData <- FindClusters(SpatialData, verbose = FALSE)
SpatialData <- RunUMAP(SpatialData, reduction = "pca", dims = 1:30)

p1 <- DimPlot(SpatialData, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(SpatialData, label = TRUE, label.size = 3)
p1 + p2

