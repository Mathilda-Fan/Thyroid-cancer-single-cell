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


