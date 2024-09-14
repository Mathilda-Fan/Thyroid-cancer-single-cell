#Identification of distinctive genes in cell clusters
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(SingleCellExperiment)
library(Matrix)
library(harmony)
library(SeuratWrappers)
library(future)
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