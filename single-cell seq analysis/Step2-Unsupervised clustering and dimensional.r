library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)

scRNA_harmony = NormalizeData(scRNA_harmony)
scRNA_harmony <- SCTransform(scRNA_harmony) %>% RunPCA() 
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "sample_id",
                            assay.use = "SCT", max.iter.harmony = 10)  

#after harmony integration
#dimension
plot2 <- ElbowPlot(scRNA_harmony, ndims=50, reduction="pca") 
pc.num=1:25
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=0.5)   ###
#umap
color_ct=c(brewer.pal(12, "Set3")[-c(2,3,9,12)],
           brewer.pal(5, "Set1")[2],
           brewer.pal(3, "Dark2")[1],
           "#fc4e2a","#fb9a99","#f781bf","#e7298a")
plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T,raster=FALSE)
plot1 = DimPlot(scRNA_harmony, reduction = "umap",split.by = 'sample_id', label=T,raster=FALSE)
plot1 = DimPlot(scRNA_harmony, reduction = "umap",group.by = 'grade', label=T,raster=FALSE)
plot1 = DimPlot(scRNA_harmony, reduction = "umap",group.by = 'Annotation', label=T,raster=FALSE)

plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = 'nFeature_RNA', 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = 'nCount_RNA', 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = 'percent.mt', 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE)
plot1 = DimPlot(scRNA_harmony, reduction = "umap",group.by = 'Phase', label=T,raster=FALSE)

#tsne
color_ct=c(brewer.pal(12, "Set3")[-c(2,3,9,12)],
           brewer.pal(5, "Set1")[2],
           brewer.pal(3, "Dark2")[1],
           "#fc4e2a","#fb9a99","#f781bf","#e7298a")
scRNA_harmony = RunTSNE(scRNA_harmony, dims =1:25, check_duplicates = FALSE)
embed_tsne <- Embeddings(scRNA_harmony, 'tsne')
write.csv(embed_tsne,'embed_tsne.csv')
plott1 = DimPlot(scRNA_harmony, reduction = "tsne") 
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "tsne", 
                      features = 'nFeature_RNA', 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "tsne", 
                      features = 'nCount_RNA', 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "tsne", 
                      features = 'percent.mt', 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE)
plot1 = DimPlot(scRNA_harmony, reduction = "tsne",group.by = 'grade', label=F,raster=FALSE)

#Tcell
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("CD3D", "CD3E",'CD8A','CD8B'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("CD4", "CD40LG",'IL7R','FOXP3'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("IL2RA", "GNLY",'NKG7','CD8A'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
#BCELL
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("MS4A1", "CD79A",'JCHAIN','PTPRC'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
#Myeloid
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("CSF1R",'FCER1G','CD68','CD163'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("S100A8",'S100A9','S100A12','FCGR3A'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("LST1",'LILRB2','CD1C','CLEC10A'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("CLEC9A",'IDO1','CXCL8','LAMP3'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
#Endothelium
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("RAMP2",'VWF','PTPRB','CLDN5'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("PECAM1",'TIMP3','ID1','KDR'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
#Fibroblasts
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("PDGFRB",'COL1A2','RGS5','TAGLN'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("COL1A1",'ACTA2','DCN','FAP'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
#EpithelialCells
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("KRT19",'CLU','TG','TPO'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
plotct1<- FeaturePlot(scRNA_harmony, 
                      reduction = "umap", 
                      features = c("FN1",'MGST1','IYD','TFF3'), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = TRUE,
                      label.size = 1)
                      
ept<-subset(scRNA_harmony, subset = celltype == "Thyrocyte")
saveRDS(ept,'Thyrocyte.rds')
fib<-subset(scRNA_harmony, subset = celltype == "Fibroblast")
saveRDS(fib,'Fibroblast.rds')
tc <-subset(scRNA_harmony, subset = celltype == "T/NK_Cell")
saveRDS(tc,'T_NK_Cell.rds')
Myeloid <-subset(scRNA_harmony, subset = celltype == "Myeloid_Cell")
saveRDS(Myeloid,'Myeloid.rds')
Endothelium <-subset(scRNA_harmony, subset = celltype == "Endothelial_Cell")
saveRDS(Endothelium,'Endothelium.rds')
Progenitor_cell <-subset(scRNA_harmony, subset = celltype == "Progenitor_Cell")
saveRDS(Progenitor_cell,'Progenitor.rds')
MAST<-subset(scRNA_harmony, subset = celltype == "Mast_Cell")
saveRDS(MAST,'Mast_Cell.rds')