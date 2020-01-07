#devtools::install_github("velocyto-team/velocyto.R")
#devtools::install_github('satijalab/seurat-wrappers')
library(Seurat)
library(SeuratWrappers)
load('~/Desktop/Monday Meeting/Data2/04skeljupyter.RData')
rm(list=ls()[-c(44)])
Sobj <- Sobj$gene

library(SingleCellExperiment)#raw counts
counts <- as.matrix(Sobj@assays$RNA@counts)#row gene col cell
sce <- SingleCellExperiment(list(counts=counts))
sce

library(scater)
mito <- grep('^mt-',rownames(sce))
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
head(colnames(colData(sce)), 10)


# libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", 
#                           log=TRUE)
# feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", 
#                           log=TRUE)
# keep <- !(libsize.drop | feature.drop )
# data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),Remaining=sum(keep))
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
# plotHighestExprs(sce, n=50) + fontsize


ave.counts <- calcAverage(sce, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
     xlab=expression(Log[10]~"average count"))

demo.keep <- ave.counts >= 1
filtered.sce <- sce[demo.keep,]
summary(demo.keep)
num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
              xlab=expression(Log[10]~"average count"))
to.keep <- num.cells > 0
sce <- sce[to.keep,]
summary(to.keep)
library(scran)
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))
# plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
#      xlab="Library size (millions)", ylab="Size factor",
#      col=c("red", "black")[sce$Oncogene], pch=16)#can put in any column here
# legend("bottomright", col=c("red", "black"), pch=16, cex=1.2,
#        legend=levels(sce$Oncogene))
sce <- normalize(sce)
var.fit <- trendVar(sce, parametric=TRUE, use.spikes=F,
                    loess.args=list(span=0.3))
var.out <- decomposeVar(sce, var.fit)
head(var.out)
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
     ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
#points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
points(var.out$mean[chosen.genes], var.out$total[chosen.genes], col="red", pch=16)
plotExpression(sce, features=rownames(var.out)[chosen.genes]) + fontsize


chosen.genesx <- rownames(sce)[order(var.out$bio, decreasing=TRUE)[1:3000]]
goal <- read.csv('~/Downloads/top3kchosengenesSCEpkg.csv')
sum(chosen.genesx %in%  goal$x)
write.csv(chosen.genesx,file='newscegenes.csv')

# library(scran)
# sce <- computeSumFactors(sce)
# summary(sizeFactors(sce))
# sce <- normalize(sce)
# fit <- trendVar(sce, parametric=TRUE,use.spikes=F)
# decomp <- decomposeVar(sce, fit)
# top.hvgs <- order(decomp$bio, decreasing=TRUE)
# head(decomp[top.hvgs,])
# plot(decomp$mean, decomp$total, xlab="Mean log-expression", ylab="Variance")
# o <- order(decomp$mean)
# lines(decomp$mean[o], decomp$tech[o], col="red", lwd=2)
# points(fit$mean, fit$var>2, col="red", pch=16)
# 
# 
# alt.fit <- trendVar(sce, use.spikes=FALSE) 
# alt.decomp <- decomposeVar(sce, alt.fit)
# 
# 
# modelGeneVar(sce)
# hvg.pbmc.var <- getTopHVGs(sce, n=1000)
# str(hvg.pbmc.var)
# 
# sum(rownames((decomp[top.hvgs,])[1:3000,]) %in% goal$x)
# 
# 
# hvg.out <- var.out[which(var.out$FDR <= 0.07 ),]
# hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
# nrow(hvg.out)
# head(hvg.out)
# 
