#Further analysis using snRNA seq data
setwd('/home/rohan/Desktop/snRNAseq/')
library(Seurat);library(dplyr);library(patchwork)
#rename features.tsv to genes.tsv after downloading
muscle_list <- list()
for (muscle in c('Quad','Sol','Tib','Mix')){
  data <- Read10X(data.dir = muscle)
  SeuratObj <- CreateSeuratObject(counts = data, project = paste0(muscle,'_snRNA'), min.cells = 3, min.features = 200)
  muscle_list[[paste0(muscle)]] <- SeuratObj
}

#Standard Preprocessing
library(ggplot2)
for (i in c(1:4)){
  muscle_list[[i]][["percent.mt"]] <- PercentageFeatureSet(muscle_list[[i]], pattern = "^mt-")
  print(VlnPlot(muscle_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  plot1 <- FeatureScatter(muscle_list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
  plot2 <- FeatureScatter(muscle_list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +NoLegend()
  print(plot1 + plot2)
}

#combine
allsamples <- merge(muscle_list[[1]], y = c(muscle_list[[2]], muscle_list[[3]],muscle_list[[4]]), add.cell.ids = names(muscle_list), project = "snRNAmuscle")
rm(muscle_list)


allsamples <- RunPCA(allsamples,features = VariableFeatures(object = allsamples))
DimHeatmap(allsamples, dims = 1:6, cells = 500, balanced = TRUE,fast=FALSE) + theme_void() + ggtitle('allsamples') 
DimHeatmap(allsamples, dims = 7:12, cells = 500, balanced = TRUE,fast=FALSE) + theme_void() + ggtitle('allsamples') 
ElbowPlot(allsamples)+ggtitle('allsamples')

#using 10 PCs for clustering
j=10
allsamples <- allsamples %>% FindNeighbors(dims = 1:j, verbose = FALSE) %>%
  FindClusters(verbose =FALSE) %>% RunUMAP(dims = 1:j)
DimPlot(allsamples, reduction = "umap")+ggtitle('allsamples')
markers <- FindAllMarkers(object = allsamples, only.pos = TRUE, min.pct = 0.5)
write.csv(markers,file = paste0('Data3/10','_fullumap_','_10pc_markers.csv'))
top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
DoHeatmap(allsamples, features = top5$gene) + NoLegend() + ggtitle('allsamples')

saveRDS(allsamples, file = "Data3/11_allsamples.rds")
rm(SeuratObj)

totest <- c('Myh4','Myh7','Myh1','Myh2','Pdgfra','Pax7','Pecam1','Pparg','Col11a1','Ptprc','Ache','Runx1','Ttn')#ptprc is cd45
length(totest)
plots=list()
for(i in 1:12){
  plots[[i]]<- FeaturePlot(allsamples,features=totest[i])+NoLegend()
}
gridExtra::grid.arrange(grobs=plots,ncol = 3)

pdf('fullUmaps_Myh.pdf')
DimPlot(allsamples, reduction = "umap")+ggtitle('allsamples')
top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
DoHeatmap(allsamples, features = top5$gene) + NoLegend() + ggtitle('allsamples')
for(i in 1:12){
  print(plots[[i]])
}
dev.off()

#selecting muscle cells
cleansamples <- subset(allsamples,idents = c(0,1,13,5,2,3,4,9))
DimPlot(cleansamples, reduction = "umap",label = TRUE)+ggtitle('allsamples')
saveRDS(cleansamples, file = "Data3/12_cleansamples.rds")
pdf('Markers.pdf')
DimPlot(cleansamples, reduction = "umap",label = TRUE)+NoLegend()+ggtitle('allsamples')
FeaturePlot(cleansamples,features='Myh4')+NoLegend()
FeaturePlot(cleansamples,features='Myh7')+NoLegend()
FeaturePlot(cleansamples,features='Myh1')+NoLegend()
FeaturePlot(cleansamples,features='Myh2')+NoLegend()
FeaturePlot(cleansamples,features='Myh8')+NoLegend()
FeaturePlot(cleansamples,features='Pecam1')+NoLegend()
FeaturePlot(cleansamples,features='Cdh5')+NoLegend()
FeaturePlot(cleansamples,features='Smtn')+NoLegend()
FeaturePlot(cleansamples,features='Lep')+NoLegend()#no leptin 
FeaturePlot(cleansamples,features='Adipoq')+NoLegend()
dev.off()

rm(allsamples)

#rename identities and findmarkers
new.cluster.ids <- c('6','7','fast 2X','fast 2A/2X','fast 2A','fast 2B','slow','8' )
names(new.cluster.ids)<-levels(cleansamples)
cleansamples <- RenameIdents(cleansamples, new.cluster.ids)
DimPlot(cleansamples, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
muscle.de.markers <- FindMarkers(cleansamples, ident.1 = 'slow', ident.2 ='fast 2X' )
write.csv(muscle.de.markers,file='Data3/15_snRNAslowvs2X.csv')

#vln plot genes of interest
extragenes <- c('Cdh4', 'Cdkl5', 'Cntn4', 'Dscam', 'Gabbr2', 'Kirrel3', 'Lingo2', 'Lrp1',
                'L1cam', 'Nrcam', 'Ntn1', 'Ntrk3', 'Ptprt', 'Ptpro', 'Robo2', 'Sdk1', 
                'Sema5a', 'Sema6d', 'Shank2', 'Sox5', 'Tnr', 'Wwox',
                'Myh4','Myh7','Myh1','Myh2','Myh8','Pecam1','Cdh5','Smtn','Lep','Adipoq')

pdf('Data3/16_snRNAseq_violinplots')
for (i in extragenes){print(VlnPlot(cleansamples,features=i,combine=T))}
dev.off()

#gene list example
df <- data.frame(gene=rownames(cleansamples@assays$RNA@data)[order(rowMeans(as.matrix(cleansamples@assays$RNA@data)),decreasing = TRUE)],
                 mean.expr=rowMeans(as.matrix(cleansamples@assays$RNA@data))[order(rowMeans(as.matrix(cleansamples@assays$RNA@data)),decreasing = TRUE)])
head(df,n = 20)

#starting from single cell analysis seurat object from part 3 as Sobj$gene
df2 <- data.frame(gene=rownames(Sobj$gene[['SCT']]@data)[order(rowMeans(as.matrix(Sobj$gene[['SCT']]@data)),decreasing = TRUE)],
                  mean.expr=rowMeans(as.matrix(Sobj$gene[['SCT']]@data))[order(rowMeans(as.matrix(Sobj$gene[['SCT']]@data)),decreasing = TRUE)])
head(df2,n = 20)

#genes expressed in 25% of cells example (also done with 96% for scrna later)
pct25N <- which(as.numeric((rowSums(data.frame(cleansamples[['RNA']]@counts>=1)))/(ncol(cleansamples[['RNA']]@counts>=1)))>=.25)
pct25N <- cleansamples[['RNA']]@counts[pct25N,]#notnormalized
pct25N <- data.frame(gene=rownames(pct25N), meanexp=rowMeans(as.matrix(pct25N)))
pct25N <- pct25N[order(pct25N$meanexp,decreasing=TRUE),]
pct25C <- which(as.numeric((rowSums(data.frame(Sobj$gene[['RNA']]@counts>=1)))/(ncol(Sobj$gene[['RNA']]@counts)))>=.25)
pct25C <- Sobj$gene[['RNA']]@counts[pct25C,]#notnormalized
pct25C <- data.frame(gene=rownames(pct25C), meanexp=rowMeans(as.matrix(pct25C)))
pct25C <- pct25C[order(pct25C$meanexp,decreasing=TRUE),]



#load in gene lists for gene annotation 
actin <- read_xlsx('genelists/Actin cytoskeleton from panther.xlsx',col_names = FALSE)[,8]
cellcyc <- read_xlsx('genelists/Cell cycle_panther.xlsx',col_names = FALSE)[,8]
mitoc <- read_xlsx('genelists/mitochondria from panther.xlsx',col_names = FALSE)[,8]
TF <- read_xlsx('genelists/mouse TFs_transcription regulator activity_from panther.xlsx',col_names = FALSE)[,8]

pct25N$category[which(rownames(pct25N) %in% actin$`...8`)] <- 'actin'
pct25N$category[which(rownames(pct25N) %in% cellcyc$`...8`)] <- 'cellcyc'
pct25N$category[which(rownames(pct25N) %in% mitoc$`...8`)] <- 'mitoc'
pct25N$category[which(rownames(pct25N) %in% TF$`...8`)] <- 'TF'
pct25C$category[which(rownames(pct25C) %in% actin$`...8`)] <- 'actin'
pct25C$category[which(rownames(pct25C) %in% cellcyc$`...8`)] <- 'cellcyc'
pct25C$category[which(rownames(pct25C) %in% mitoc$`...8`)] <- 'mitoc'
pct25C$category[which(rownames(pct25C) %in% TF$`...8`)] <- 'TF'

dfn <- data.frame(group=names(table(pct25N$category)),value=as.numeric(table(pct25N$category)))
dfn$value <- 100*dfn$value/sum(dfn$value)
dfc <- data.frame(group=names(table(pct25C$category)),value=as.numeric(table(pct25C$category)))
dfc$value <- 100*dfc$value/sum(dfc$value)
print('snrna');dfn;print('scrna');dfc

#comparing fast 2x vs slow/fast2a t test
new.cluster.ids <- c("fast2xb" , "fast2xa", "slow/fast2A") 
names(new.cluster.ids) <- levels(Sobj$gene)
Sobj$gene <- RenameIdents(Sobj$gene, new.cluster.ids)
originalscde <- FindMarkers(Sobj$gene,ident.1 = 'slow/fast2A',ident.2 = 'fast2x')
originalscde <- FindMarkers(Sobj$gene,ident.1 = 'slow/fast2A',ident.2 = 'fast2x',test.use = 't')
write.csv(originalscde, file='Data3/28_slowplusfast2a_vs_2x_TTESTinSCRNA.csv')

shortsn <- subset(cleansamples,idents = c('slow','fast 2A','fast 2X'))
new.cluster.ids <- c("fast2x","slow/fast2a" , "slow/fast2a") 
names(new.cluster.ids) <- levels(shortsn)
shortsn <- RenameIdents(shortsn, new.cluster.ids)
DimPlot(shortsn, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(Sobj$gene, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

snde <- FindMarkers(shortsn,ident.1 = 'slow/fast2a',ident.2 = 'fast2x')
write.csv(snde, file='Data3/25_slowplusfast2a_vs_2x_inSNRNA.csv')