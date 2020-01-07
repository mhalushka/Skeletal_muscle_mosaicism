
#load Libraries
library(Seurat);library(dplyr);library(ggplot2)
load(file="/home-4/rverma6@jhu.edu/muscle_proj/Data/009skeljupyter.RData")
ls()

#read in and use sctransform to clean
gene <- CreateSeuratObject(counts = pbmc@raw.data, project = "skel_musc", min.cells = 3, min.features = 200)
gene[["percent.mt"]] <- PercentageFeatureSet(gene, pattern = "^mt-")
gene@meta.data[['myocyte_type']] <- pbmc@meta.data$celltype
gene <- subset(gene, subset = nFeature_RNA > 5000 & percent.mt < 20)
gene <- SCTransform(gene, vars.to.regress = "percent.mt", verbose = FALSE)

library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = manngenes$Gene.names , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
manngenesmouse <- (genesV2[, 2])
length(unique(manngenes$Gene.names));length(manngenesmouse)
manngenes <- gdata::read.xls('/home-4/rverma6@jhu.edu/muscle_proj/Data/Mann Cell Reports data.xlsx',sheet=1,skip=1)
manngenes$Ratio.of.slow.to.fast[-c(1:4)] <- log(manngenes$Ratio.of.slow.to.fast[-c(1:4)])#log transform all those values
manngeneadj <- manngenes$Ratio.of.slow.to.fast[-c(1:4)]
cutoffgenes <- c(c(1:4),which(manngeneadj<mean(manngeneadj)-1*sd(manngeneadj)),which(manngeneadj>mean(manngeneadj)+1*sd(manngeneadj)))
manngenesmouse <- manngenes$Gene.names[cutoffgenes]
length(manngenesmouse)
manngenesmouse
genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = manngenesmouse , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T) 
manngenesmouse <- unique(genesV2[, 2])
manngenesmouse

#New objects for Dimensional Reductions; PCA
protein <- gene ; mann <- gene
gene <- RunPCA(gene, features = VariableFeatures(object = gene))
protein <- RunPCA(protein, features = c(mosaiclistfiltered$Mouse.Gene))#recall filter removed unknown and NA q vals
mann <- RunPCA(mann, features = manngenesmouse)

Sobj <- list(gene,protein,mann)
names(Sobj) <- c('gene','protein','mann')

for(i in 1:length(Sobj)){
    write.csv(Loadings(Sobj[[i]], reduction = "pca"),file=paste0('/home-net/home-4/rverma6@jhu.edu/muscle_proj/Data2/0',i,'_',names(Sobj)[i],'PCdata.csv'))
    print(DimHeatmap(Sobj[[i]], dims = 1:6, cells = 500, balanced = TRUE,fast=FALSE) + theme_void() + ggtitle(names(Sobj)[i])) 
    print(DimHeatmap(Sobj[[i]], dims = 7:12, cells = 500, balanced = TRUE,fast=FALSE) + theme_void() + ggtitle(names(Sobj)[i])) 
    print(ElbowPlot(Sobj[[i]])+ggtitle(names(Sobj)[i]))
}

#saw doing 20pcs generated good clustering and sctransform effectively helps out want to see 
head(Stdev(mann, reduction = "pca"),20)#20 is at 1.15 and 32 is just above 1

head(Stdev(gene, reduction = "pca"),20)#20 is at 1.15 and 32 is just above 1

head(Stdev(protein, reduction = "pca"),20)#20 is at 1.15 and 32 is just above 1

for (i in 1:length(Sobj)){
    if(i==1){j<-4};if(i==2){j<-4};if(i==3){j<-4}#set pc's
    tmp <- Sobj[[i]] %>% RunUMAP(dims = 1:j) %>% FindNeighbors(dims = 1:5, verbose = FALSE) %>%
        FindClusters(verbose =FALSE)
    print(DimPlot(tmp, reduction = "umap")+ggtitle(names(Sobj)[i]))
    markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.5)
    write.csv(markers,file = paste0('/home-net/home-4/rverma6@jhu.edu/muscle_proj/Data2/0',i+4,'umap_',names(Sobj)[i],'_05pc_markers.csv'))
    top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
    print(DoHeatmap(tmp, features = top5$gene) + NoLegend() + ggtitle(names(Sobj)[i]))
    Sobj[i]<- tmp
}

Sobj

Sobj[[1]]

#add metadata clusters
Sobj[[1]]@meta.data$gene <- Sobj[[1]]@meta.data$SCT_snn_res.0.8
Sobj[[1]]@meta.data$protein <- Sobj[[2]]@meta.data$SCT_snn_res.0.8
Sobj[[1]]@meta.data$mann <- Sobj[[3]]@meta.data$SCT_snn_res.0.8

Sobj[[2]]@meta.data$gene <- Sobj[[1]]@meta.data$SCT_snn_res.0.8
Sobj[[2]]@meta.data$protein <- Sobj[[2]]@meta.data$SCT_snn_res.0.8
Sobj[[2]]@meta.data$mann <- Sobj[[3]]@meta.data$SCT_snn_res.0.8

Sobj[[3]]@meta.data$gene <- Sobj[[1]]@meta.data$SCT_snn_res.0.8
Sobj[[3]]@meta.data$protein <- Sobj[[2]]@meta.data$SCT_snn_res.0.8
Sobj[[3]]@meta.data$mann <- Sobj[[3]]@meta.data$SCT_snn_res.0.8

gridExtra::grid.arrange(
    DimPlot(Sobj$gene, reduction = "umap",label = F,group.by='gene')+ggtitle('Gene')+NoLegend(),
    DimPlot(Sobj$protein, reduction = "umap",label = F,group.by='gene')+ggtitle('Protein')+NoLegend(),
    DimPlot(Sobj$mann, reduction = "umap",label = F,group.by='gene')+ggtitle('Mann')+NoLegend(),
    
    DimPlot(Sobj$gene, reduction = "umap",label = F,group.by='protein')+NoLegend(),
    DimPlot(Sobj$protein, reduction = "umap",label = F,group.by='protein')+NoLegend(),
    DimPlot(Sobj$mann, reduction = "umap",label = F,group.by='protein')+NoLegend(),
    
    DimPlot(Sobj$gene, reduction = "umap",label = F,group.by='mann')+NoLegend(),
    DimPlot(Sobj$protein, reduction = "umap",label = F,group.by='mann')+NoLegend(),
    DimPlot(Sobj$mann, reduction = "umap",label = F,group.by='mann')+NoLegend(),
    nrow = 3,
    top = "UMAP Comparisons",
    bottom = grid::textGrob(
    "Row 1 colored by gene 2 by protein 3 by mann",
    gp = grid::gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  )
)

geness <- VariableFeatures(gene);length(geness)
geness
protss <- (c(mosaiclistfiltered$Mouse.Gene));length(protss)
protss
mannss <- manngenesmouse;length(mannss)
mannss
library(VennDiagram)
venn.diagram(list('gene'=geness,'prot'=protss,'mann'=mannss),fill = c("yellow","cyan","red"), cex = 1.5, filename="/home-4/rverma6@jhu.edu/muscle_proj/Data2/08redo2venndiagram_all3.tiff")

intersect(c(intersect(mannss,geness)),protss)
intersect(mannss,protss)
intersect(mannss,geness)
intersect(geness,protss)
setdiff(mannss,protss)#all in mannss not in protss
setdiff(mannss,c(protss,geness))#all in mannss not in protss/genes
setdiff(geness,c(protss,mannss))#genes not in either
setdiff(protss,c(geness,mannss))

save(list=ls(),file='/home-4/rverma6@jhu.edu/muscle_proj/Data2/04skeljupyter.RData')

load('/home-4/rverma6@jhu.edu/muscle_proj/Data2/04skeljupyter.RData')

library(Seurat);library(ggplot2)

colors <- scales::hue_pal()(12)
colors

gridExtra::grid.arrange(
    DimPlot(Sobj$gene, reduction = "umap",label = F,group.by='gene',cols = colors[1:3])+ggtitle('Gene')+NoLegend(),
    DimPlot(Sobj$protein, reduction = "umap",label = F,group.by='gene',cols = colors[1:3])+ggtitle('Protein')+NoLegend(),
    DimPlot(Sobj$mann, reduction = "umap",label = F,group.by='gene',cols = colors[1:3])+ggtitle('Mann')+NoLegend(),
    
    DimPlot(Sobj$gene, reduction = "umap",label = F,group.by='protein',cols = colors[4:7])+NoLegend(),
    DimPlot(Sobj$protein, reduction = "umap",label = F,group.by='protein',cols = colors[4:7])+NoLegend(),
    DimPlot(Sobj$mann, reduction = "umap",label = F,group.by='protein',cols = colors[4:7])+NoLegend(),
    
    DimPlot(Sobj$gene, reduction = "umap",label = F,group.by='mann',cols = colors[8:12])+NoLegend(),
    DimPlot(Sobj$protein, reduction = "umap",label = F,group.by='mann',cols = colors[8:12])+NoLegend(),
    DimPlot(Sobj$mann, reduction = "umap",label = F,group.by='mann',cols = colors[8:12])+NoLegend(),
    nrow = 3,
    top = "UMAP Comparisons",
    bottom = grid::textGrob(
    "Row 1 colored by gene 2 by protein 3 by mann",
    gp = grid::gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  )
)

pdf('../Data2/10umapcomparisions.pdf')
gridExtra::grid.arrange(
    DimPlot(Sobj$gene, reduction = "umap",label = F,group.by='gene',cols = colors[1:3])+ggtitle('Gene')+NoLegend(),
    DimPlot(Sobj$protein, reduction = "umap",label = F,group.by='gene',cols = colors[1:3])+ggtitle('Protein')+NoLegend(),
    DimPlot(Sobj$mann, reduction = "umap",label = F,group.by='gene',cols = colors[1:3])+ggtitle('Mann')+NoLegend(),
    
    DimPlot(Sobj$gene, reduction = "umap",label = F,group.by='protein',cols = colors[4:7])+NoLegend(),
    DimPlot(Sobj$protein, reduction = "umap",label = F,group.by='protein',cols = colors[4:7])+NoLegend(),
    DimPlot(Sobj$mann, reduction = "umap",label = F,group.by='protein',cols = colors[4:7])+NoLegend(),
    
    DimPlot(Sobj$gene, reduction = "umap",label = F,group.by='mann',cols = colors[8:12])+NoLegend(),
    DimPlot(Sobj$protein, reduction = "umap",label = F,group.by='mann',cols = colors[8:12])+NoLegend(),
    DimPlot(Sobj$mann, reduction = "umap",label = F,group.by='mann',cols = colors[8:12])+NoLegend(),
    nrow = 3,
    top = "UMAP Comparisons",
    bottom = grid::textGrob(
    "Row 1 colored by gene 2 by protein 3 by mann",
    gp = grid::gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  )
)
dev.off()

FeaturePlot(Sobj$gene,features='Actn3')
FeaturePlot(Sobj$protein,features='Actn3')
FeaturePlot(Sobj$mann,features='Actn3')

VlnPlot(Sobj$gene, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = 'orig.ident' )

ls()
