#Pathway enrichment analysis
#GSVA
library(msigdbr)
library(org.Hs.eg.db)
#genesets <- msigdbr(species = 'Homo sapiens',category = 'C2',subcategory = 'CP:KEGG')
genesets <- msigdbr(species = 'Homo sapiens',category = 'C5',subcategory = 'GO:BP')
genesets <- subset(genesets,select = c('gs_name','gene_symbol'))%>%as.data.frame()
genesets <- split(genesets$gene_symbol,genesets$gs_name)
expr <- AverageExpression(Epithelium,assays = 'RNA',layer = 'data')[[1]]
expr <- expr[rowSums(expr)>0.1,]
expr <- as.matrix(expr)
head(expr)
gsva.res <- gsva(expr,genesets,method='gsva')
GSVA.ID <- str_replace(row.names(gsva.res),'GO_','')
gsva.res <- data.frame(ID=GSVA.ID,gsva.res)
row.names(gsva.res) <- gsva.res$ID
gsva.res <- gsva.res[-1]
gsva.df <- data.frame(Genesets=rownames(gsva.res),gsva.res,check.names = F)
plot1 <- pheatmap::pheatmap(gsva.res[1701:1713,],show_colnames = T,scale = 'row',angle_col = '45',fontsize_number = 0.1,color = colorRampPalette(c('navy','white','firebrick3'))(50))


