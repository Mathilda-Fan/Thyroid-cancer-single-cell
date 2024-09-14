#SCENIC Analysis
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(harmony)

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
subcell.id <- sample(colnames(CAF),1500)
scRNAsub <- scRNA_harmony[,subcell.id]
saveRDS(scRNAsub, "scRNAsub.rds")
scRNAsub=readRDS("scRNAsub.rds")

exprMat <- GetAssayData(object = scRNAsub, layer = "counts")
exprMat <- as.matrix(exprMat1)
mydbDIR <- "D:\\single-wzynew\\TF"
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=14,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "os")
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
interestingGenes <- c("PROS1")
interestingGenes[which(!interestingGenes %in% genesKept)]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)

runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
head("int/1.4_GENIE3_linkList.Rds")
head(readRDS("int/1.4_GENIE3_linkList.Rds") )

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) 
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 

scRNA_harmony[["RNA3"]] <- as(object = scRNA_harmony[["RNA"]], Class = "Assay")
DefaultAssay(scRNA_harmony) <- "RNA3"
scRNA_harmony[["RNA"]] <- NULL
scRNA_harmony <- RenameAssays(object = scRNA_harmony, RNA3 = 'RNA')

scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "os")
library(foreach)
exprMat_all <- as.matrix(scRNA_harmony@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
rm(scRNA_harmony)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)
?runSCENIC_3_scoreCells

runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all) 

##meta
cellInfo <- data.frame(scRNA_harmony@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="Annotation1")] <- "Annotation"
cellInfo <- cellInfo[,c("sample","cluster","Annotation1")]
saveRDS(cellInfo, file="int/cellInfo.Rds")

mydbDIR <- "D:\\single-wzynew\\TF"
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "os")


library(foreach)
nPcs <- c(5)
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
par(mfrow=c(length(nPcs), 1))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="Annotation1",   cex=.5)
?plotTsne_compareSettings
?AUCell::plotTsne_cellProps
par(mfrow=c(2,2))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="Annotation1208", cex=.5)
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 5
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(scRNA_harmony, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')

BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA_harmony, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'scRNAbin.rds')

scRNA_harmony = RunTSNE(scRNA_harmony, dims =1:25, check_duplicates = FALSE)
embed_tsne <- Embeddings(scRNA_harmony, 'tsne')
plott1 = DimPlot(scRNA_harmony, reduction = "tsne") 
plott1

dir.create('scenic_seurat')

library(pheatmap)
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select = 'Annotation1208')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
#regulons
rownames(AUCmatrix)[1:1000]
regulonsNames <- c("FOXP2_extended_437g","FOXP2_228g","NFYB_extended_1290g","NFYB_1051g","ETS1_extended_6164g","ETS1_3787g","KDM5A_736g","KDM5A_extended_4779g","PML_extended_1756g","PML_947g")
my.regulons <- regulonsNames
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]
plot1<- pheatmap(myAUCmatrix, show_colnames=F, annotation_col=celltype )
plot1 <- pheatmap(myBINmatrix, show_colnames=F, annotation_col=celltype,
         color = colorRampPalette(colors = c("white","black"))(100))

