#Quality control and batch effect mitigation in scRNA-Seq data analysis
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(SingleCellExperiment)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(Rcpp)
library(harmony)
library(SeuratData)
library(celldex)
library(scRNAseq)
library(data.table)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(clusterProfiler)
library(SeuratWrappers)
library(stringr)
library(RColorBrewer)
library(scales)
library(reshape2)
library(future)

setwd('D:\\single-wzynew\\data')
files=dir("D:\\single-wzynew\\data")
id <-dir("D:\\single-wzynew\\data")
scRNAlist <- list()
for(i in 1:length(files)){
  count <-Read10X(data.dir = files[i])
  scRNAlist[[i]] = CreateSeuratObject(counts = count,min.cells = 3,min.features = 200)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]][["sample_id"]] <-id[i]
  scRNAlist[[i]] = RenameCells(scRNAlist[[i]], add.cell.id = id[i])
}
scRNA_harmony<-merge(x=scRNAlist[[1]],y=c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],
                                          scRNAlist[[5]],scRNAlist[[6]],scRNAlist[[7]],scRNAlist[[8]],scRNAlist[[9]],scRNAlist[[10]],scRNAlist[[11]],scRNAlist[[12]],scRNAlist[[13]],scRNAlist[[14]]))

GroupA<-c('A1','A2','A3','A4')
GroupB<-c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10')
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$sample_id %in% GroupA),'grade'] <- 'GroupA'
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$sample_id %in% GroupB),'grade'] <- 'GroupB'
lev <- c('A1','A2','A3','A4','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10')
scRNA_harmony$sample_id=factor(scRNA_harmony$sample_id,levels = lev)
mito_plot <- VlnPlot(scRNA_harmony,
                     features ="percent.mt",
                     group.by ="sample_id",
                     pt.size = 0.001, 
                     ncol = 1) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
nFeature_plot <- VlnPlot(scRNA_harmony,
                         features ="nFeature_RNA",
                         group.by ="sample_id",
                         pt.size = 0.01, 
                         ncol = 1) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
nCount <- VlnPlot(scRNA_harmony,
                  features ="nCount_RNA",
                  group.by ="sample_id",
                  pt.size = 0.01, 
                  ncol = 1) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
plot11=FeatureScatter(scRNA_harmony, feature1 = "nCount_RNA", feature2 = "percent.mt")+ scale_y_continuous(limits = c(0,100))+scale_x_continuous(limits = c(0,1000000))
plot21=FeatureScatter(scRNA_harmony, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ scale_y_continuous(limits = c(0,15000))+scale_x_continuous(limits = c(0,1000000))
pearplot1 <- CombinePlots(plots = list(plot11, plot21), nrow=1, legend="none")

# PercentageFeatureSet
scRNA_harmony <- PercentageFeatureSet(scRNA_harmony, "^RP[SL]", col.name = "percent_ribo")
scRNA_harmony <- PercentageFeatureSet(scRNA_harmony, "^HB[^(P)]", col.name = "percent_hb")
feats <- c("percent_ribo", "percent_hb")
# VlnPlot
VlnPlot(scRNA_harmony, group.by = "sample_id", features = feats, pt.size = 0.1, ncol = 3) + 
  NoLegend()
ected_f <- rownames(scRNA_harmony)[Matrix::rowSums(scRNA_harmony) > 3]
scRNA_harmony <- subset(scRNA_harmony, features = selected_f)
#filter
scRNA_harmony <- subset(scRNA_harmony, subset = nFeature_RNA > 200& nFeature_RNA<6000 &percent.mt < 10 & nCount_RNA >2000&nCount_RNA <40000) 
dim(scRNA_harmony)

# Compute the relative expression of each gene per cell Use sparse matrix
# operations, if your dataset is large, doing matrix devisions the regular way
# will take a very long time.
par(mar = c(4, 8, 2, 1))
C <- scRNA_harmony@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell", 
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

dim(scRNA_harmony)
## [1] 18147  5762
# Filter Mitocondrial
scRNA_harmony <- scRNA_harmony[!grepl("^MT-", rownames(scRNA_harmony)), ]
# Filter Ribossomal gene (optional if that is a problem on your data) scRNA_harmony
scRNA_harmony <- scRNA_harmony[ ! grepl('^RP[SL]', rownames(scRNA_harmony)), ]
# Filter Hemoglobin gene (optional if that is a problem on your data)
scRNA_harmony <- scRNA_harmony[!grepl("^HB[^(P)]", rownames(scRNA_harmony)), ]
dim(scRNA_harmony)
## [1] 18121  5762
# Before running CellCycleScoring the data need to be normalized and logtransformed.
scRNA_harmony = NormalizeData(scRNA_harmony)

# CellCycleScoring
scRNA_harmony <- CellCycleScoring(object = scRNA_harmony, g2m.features = cc.genes$g2m.genes, 
                                  s.features = cc.genes$s.genes)
VlnPlot(scRNA_harmony, features = c("S.Score", "G2M.Score"), group.by = "sample_id", 
        ncol = 4, pt.size = 0.1)
#number of cell after filtration
cell.num= as.data.frame(table(scRNA_harmony$sample_id))
colnames(cell.num)<-c('Sample id','Number of cells')

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

#FindAllMarkers-Top30
diff.wilcox = FindAllMarkers(scRNA_harmony)
all.markers = diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05,)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)
write.csv(top30, "top30_diff_genes_wilcox.csv", row.names = F)
marker_celltype=FindAllMarkers(scRNA_harmony,only.pos = T)
library(xlsx)
marker_celltype=marker_celltype%>%filter(p_val_adj < 0.05)
marker_celltype$d=marker_celltype$pct.1-marker_celltype$pct.2
marker_celltype=marker_celltype%>%filter(d > 0.3)
marker_celltype=marker_celltype%>%arrange(cluster,desc(avg_log2FC))
marker_celltype=as.data.frame(marker_celltype)
write.xlsx(marker_celltype,file = "markers_scRNA_harmony.xlsx",row.names = F)

 select_genes<-c('EPCAM','CLU','TG',
                'CD3E','IL7R','CD3D',
                'MS4A1','CD79A','IGKC',
                'FCER1G','CSF1R','LST1',
                'COL1A2','PDGFRB','ACTA2',
                'RAMP2','VWF','PTPRB')
p <-DotPlot(scRNA_harmony, features = select_genes ) 

scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters %in% c('0','6','14','16','17','19','21')),'Annotation'] <- 'Thyrocyte'
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters %in% c('5','11')),'Annotation'] <- 'Myeloid_Cell'
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters %in% c('3','12','20')),'Annotation'] <- 'Endothelial_Cell'
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters %in% c('2','7','9','10')),'Annotation'] <- 'T/NK_Cell'
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters %in% c('4','18')),'Annotation'] <- 'B_Cell'
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters %in% c('13')),'Annotation'] <- 'Plasma_Cell'
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters %in% c('1','8','22')),'Annotation'] <- 'Fibroblast'
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters %in% c('15')),'Annotation'] <- 'Progenitor_Cell'
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters %in% c('23')),'Annotation'] <- 'Mast_Cell'
plot2 = DimPlot(scRNA_harmony, reduction = "umap",group.by = 'Annotation', label=F,raster=FALSE,cols = color_ct)
plot2 = DimPlot(scRNA_harmony, reduction = "tsne",group.by = 'Annotation', label=F,raster=FALSE,cols = color_ct)

scRNA_harmony[['celltype']]<-scRNA_harmony$Annotation

##cell number by celltype
prop.table(table(scRNA_harmony$Annotation))
cell.prop<-as.data.frame(prop.table(table(scRNA_harmony$Annotation, scRNA_harmony$sample_id)))
colnames(cell.prop)<-c('celltype',"sample_id","proportion")
p<-ggplot(cell.prop,aes(sample_id,proportion,fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))
