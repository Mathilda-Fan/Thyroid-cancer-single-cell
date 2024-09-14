#Cellular communication analysis
library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)
library(xlsx)
library(harmony)

load("Thyrocyte.Rdata")
load("mianyixibao.Rdata")
load("Fibroblast.Rdata")
load("Endothelial_cell.Rdata")
Thyrocyte$cluster_id <-paste('C',Thyrocyte$seurat_clusters,sep='') 
Thyrocyte_C5 <- subset(Thyrocyte, subset = cluster_id == "C5")
saveRDS(Thyrocyte_C5,'Thyrocyte_C5.rds')
mianyixibao$cluster_id <- mianyixibao$AnnotationNEW
Fibroblast$cluster_id <- Fibroblast$AnnotationNEW
Endothelial_cell$cluster_id <- Endothelial_cell$AnnotationNEW
Idents(Thyrocyte_C5) <-Thyrocyte_C5$cluster_id
Idents(mianyixibao) <-mianyixibao$cluster_id
Idents(Fibroblast) <-Fibroblast$cluster_id
Idents(Endothelial_cell) <-Endothelial_cell$cluster_id

scrna<-merge(x=Thyrocyte_C5,y=c(mianyixibao,Fibroblast,Endothelial_cell))
saveRDS(scrna,'Thyrocyte_C5_and_othertype_cellchat.rds')
scrna <- `merge_thy_DC_cellchat——grade`
rm(`merge_thy_DC_cellchat——grade`)
###cellchat
table(scrna$cluster_id)
data.input <- GetAssayData(scrna, assay = "RNA", slot = "data")
identity <- subset(scrna@meta.data, select = "cluster_id")
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
table(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact",key = "annotation")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 15)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,type = "truncatedMean",trim = 0.1)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat@netP$pathways
df.net <- subsetCommunication(cellchat)
write.csv(df.net,file = 'cellchat-Cell-Cell.csv')
table(cellchat@idents)
df.net <- subsetCommunication(cellchat, sources.use = c(4:5), targets.use = c(1:3))
groupSize <- as.numeric(table(cellchat@idents))
pdf('interaction_Thyrocyte_C5_and_othertype_Cell-Cell.pdf')
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
pdf('individual_Thyrocyte_C5_and_othertype_Cell-Cell.pdf')
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
p<-netVisual_bubble(cellchat, sources.use = c(1:7), targets.use = c(11), remove.isolate = FALSE)
p1<-netVisual_bubble(cellchat, sources.use = c(4:5), targets.use = c(1:3), remove.isolate = FALSE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pathways.show <- c("CXCL") 
netVisual_aggregate(cellchat, signaling = pathways.show)
#heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, color.heatmap = "Reds")
pathways.show.all <- cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)


