
#load Libraries
library(Seurat);library(dplyr);library(ggplot2)
load(file="/home-4/rverma6@jhu.edu/muscle_proj/Data2/04skeljupyter.RData")
ls()

gene <- Sobj[[1]];protein <- Sobj[[2]]; mann <- Sobj[[3]]

#checking fragment locations
DimPlot(gene, reduction = "umap",label = F,group.by='myocyte_type')+ggtitle('Gene')+NoLegend()
DimPlot(protein, reduction = "umap",label = F,group.by='myocyte_type')+ggtitle('Protein')+NoLegend()
DimPlot(mann, reduction = "umap",label = F,group.by='myocyte_type')+ggtitle('Mann')+NoLegend()

genestotest=c( 'Ttn','Xirp2','Dmd','Fhl1', 'Des', 'Zbtb20', 'Tmod4')
for (i in genestotest){print(FeaturePlot(object = gene, features = i))}
for (i in genestotest){print(VlnPlot(object = gene, features = i))}

for (i in genestotest){print(FeaturePlot(object = protein, features = i))}
for (i in genestotest){print(VlnPlot(object = protein, features = i))}

for (i in genestotest){print(FeaturePlot(object = mann, features = i))}
for (i in genestotest){print(VlnPlot(object = mann, features = i))}

genestotestdf <- gene@assays$SCT@data[genestotest,]
genestotestdf <- t(as.data.frame(genestotestdf))
genestotestdf <- as.data.frame(genestotestdf)
for (i in genestotest){print(ggplot(data = genestotestdf,aes(x=get(i)))+geom_density()+ggtitle(i))}

#match(rownames(Sobj[[i]]@meta.data),rownames(pbmc@meta.data))
#rownames(Sobj[[1]]@meta.data)[9]==rownames(pbmc@meta.data)[21]
for(i in 1:3){
    total_reads_i <- match(rownames(Sobj[[i]]@meta.data),rownames(pbmc@meta.data))
    Sobj[[i]]@meta.data$total_reads <- pbmc@meta.data$total_reads[total_reads_i]
}

for (i in 1:length(Sobj)){
    dfe<- as.data.frame(Embeddings(Sobj[[i]][["umap"]]))
    print(ggplot(dfe,aes(x=UMAP_1,y=UMAP_2))+
    geom_point(aes(color = (Sobj[[i]]@meta.data$total_reads)))+
    scale_colour_gradientn(colours=rainbow(4))+ggtitle(names(Sobj)[i]))
    print(ggplot(dfe,aes(x=UMAP_1,y=UMAP_2))+geom_point(aes(color = Sobj[[i]]@meta.data$percent.mt))+
    scale_colour_gradientn(colours=rainbow(4))+ggtitle(names(Sobj)[i]))
}

for (i in 1:3) {print(VlnPlot(Sobj[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))}

#extra images gene, protein, mann
manngenes.x <- gdata::read.xls('/home-4/rverma6@jhu.edu/muscle_proj/Data/Mann Cell Reports data.xlsx',sheet=1,skip=1)
which(manngenes$Ratio.of.slow.to.fast.1==0)
range(manngenes$Ratio.of.slow.to.fast.1==+Inf)

library(biomaRt)
manngenes$Gene.names[which(manngenes.x$Ratio.of.slow.to.fast==0|manngenes.x$Ratio.of.slow.to.fast>100)]
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = manngenes$Gene.names[which(manngenes.x$Ratio.of.slow.to.fast==0|manngenes.x$Ratio.of.slow.to.fast>100)] , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
(genesV2[, 2])

for (i in genesV2[, 2]){print(FeaturePlot(object = gene, features = i))}
for (i in genesV2[, 2]){print(VlnPlot(object = gene, features = i))}

for (i in genesV2[, 2]){print(FeaturePlot(object = protein, features = i))}
for (i in genesV2[, 2]){print(VlnPlot(object = protein, features = i))}

for (i in genesV2[, 2]){print(FeaturePlot(object = mann, features = i))}
for (i in genesV2[, 2]){print(VlnPlot(object = mann, features = i))}

intersecting <- intersect(c(intersect(mannss,geness)),protss)
for (i in intersecting){print(FeaturePlot(object = gene, features = i))}
for (i in intersecting){print(VlnPlot(object = gene, features = i))}

for (i in intersecting){print(FeaturePlot(object = protein, features = i))}
for (i in intersecting){print(VlnPlot(object = protein, features = i))}

for (i in intersecting){print(FeaturePlot(object = mann, features = i))}
for (i in intersecting){print(VlnPlot(object = mann, features = i))}

#Requested Histograms for genes not used by gene umap clustering
questionable1 <- c('Aldoa', 'Eno3', 'Gapdh', 'Ldha', 'Mybpc2', 'Mylpf', 'Pfkm', 'Pkm', 'Tnni2', 'Tnnt3', 'Tpm1')

questionable <- gene@assays$SCT[questionable1,]
questionable <- t(as.data.frame(questionable))
questionable <- as.data.frame(questionable)
for (i in questionable1){print(ggplot(data = questionable,aes(x=get(i)))+geom_density()+ggtitle(i))}

forhist1 <- intersect(protss,mannss)[!intersect(protss,mannss) %in% intersect(protss,intersect(mannss,geness))]
forhist <- gene@assays$SCT[forhist1,]
forhist <- t(as.data.frame(forhist))
forhist <- as.data.frame(forhist)
for (i in forhist1){print(ggplot(data = forhist,aes(x=get(i)))+geom_density()+ggtitle(i))}

addngenes1 <- c('Myh1', 'Myh8', 'Myl3' , 'Myoz2')
addngenes <- gene@assays$SCT[addngenes1,]
addngenes <- t(as.data.frame(addngenes))
addngenes <- as.data.frame(addngenes)
for (i in addngenes1){print(ggplot(data = addngenes,aes(x=get(i)))+geom_density()+ggtitle(i))}

pdf('/home-4/rverma6@jhu.edu/muscle_proj/Data2/09_Myh1_8_Myl3_Myoz2_histograms.pdf')
for (i in addngenes1){print(ggplot(data = addngenes,aes(x=get(i)))+geom_density()+ggtitle(i))}
dev.off()
