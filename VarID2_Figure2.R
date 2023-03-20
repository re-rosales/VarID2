# This script contains the code to reproduce the analysis of nuclear and cellular noise transcripts with VarID2
# Contained in Figures 2 and S2 of Rosales-Alvarez et al., manuscript
# Data from 10X genomics, samples from peripheral blood mononuclear cells (PBMC) from a healthy donor - granulocytes removed through cell sorting (10k)
# The code in this script was performed with RaceID v0.2.6. For further versions, cluster numbers and parameter settings potentially change.

# Setup python environment in order to use leiden clustering
# It requires to install reticulate R package and the python "leidenalg" and "igraph" modules in the indicated python environment
reticulate::use_python("./anaconda3/bin/python", required = TRUE)
# Corroborate "leidenalg" and "igraph" are available, as well as python
reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()

# Load required libraries
library(RaceID)
library(Matrix)
library(pheatmap)
library(ggplot2)
library(scales)
library(Seurat)
library(magrittr)


#Import multiomics data, whose expression matrix contains the nuclei data
#Download input data from: 10x Genomics,
#https: //support.10xgenomics.com/single-cell-multiome-atac-gex/ datasets/1.0.0/pbmc_granulocyte_sorted_10k
inputdata.10x <- Read10X_h5("/path_to_data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
d <- inputdata.10x$`Gene Expression`

#Import cell data
#Download input data from: 10x Genomics,
#https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3
d2 <- Seurat::Read10X_h5("/data/gruen/group/rosales/ATACseq/pbmc_10k_v3_filtered_feature_bc_matrix.h5")

#Comparison of UMI counts and number of detected features per cell barcode
quanL <-list(NucUMIs = colSums(d), CellUMIs = colSums(d2),
             NucGenes = apply(d,2, function(x) sum(x != 0)),
             CellGenes = apply(d2,2, function(x) sum(x != 0)) )
boxplot(quanL[1:2], main="UMIs per cell", log="y", pch=20, cex=0.2)
boxplot(quanL[3:4], main="Detected genes per cell", pch=20, cex=0.2)

#Create a SCseq object, perform quality filtering, knn inference and pruning, clustering and noise estimation
#Start with nuclei data
#Filter out barcodes with more than 25000 transcripts
f <- colSums(d) < 25000
sN <- SCseq(d[,f])
sN <- filterdata(sN,mintotal=1000,FGenes=grep("^(MT-|RP(L|S)|GM\\d)",rownames(d),value=TRUE))
expData  <- as.matrix(getExpData(sN))
res      <- pruneKnn(expData,no_cores=5)
cl    <- graphCluster(res,pvalue=0.01,leiden.resolution=2, use.leiden=TRUE)

#Noise estimation and related quantities:
x <- getFilteredCounts(sN)
noise <- compTBNoise(res,x,pvalue=0.01,gamma=1,no_cores=5)
qn <- quantKnn(res, noise, sN, pvalue = 0.01, minN = 5, no_cores = 5)

sN <- updateSC(sN,res=res,cl=cl,noise=noise,flo=.1)
sN <- updateSC(sN,res=res,cl=cl)
sN <- comptsne(sN)
sN <- compumap(sN)

#UMAP plot
plotmap(sN, um=TRUE)
#Plot of marker genes
genes <- c("CD8B","CD8A", "LEF1", "IL7R", "CCR7", "S100A4", "INPP4B","GNLY", "NKG7", "MS4A1", "IGHM", "CD14", "LYZ",
       "VCAN", "FCGR3A","MS4A7","FCER1A", "CST3", "RHEX","HLA-DPA1","CD74")
fractDotPlot(sN, genes, cluster=c(16,2,17,3,11,18,1,14,20,6,12,23,10,13,22,4,7,5,9,15,8,19,21), logscale=TRUE)
#Plot noise levels per cluster
clus <- c(16,2,17,3,11,18,1,14,20,6,12,23,10,13,22,4,7,5,9,15,8,19,21)
plotB(qn$noise.av,sN@cpart, cluster=StemCluster, set=clus, pch = 20, cex = 0.2,col = sN@fcol[clus], ylab="biological noise")


#Perform the same analysis for cell data:
#Filter out barcodes with more than 25000 transcripts
f <- colSums(d2) < 25000
sN2 <- SCseq(d2[,f])
sN2 <- filterdata(sN2,mintotal=1000,FGenes=grep("^(MT-|RP(L|S)|GM\\d)",rownames(d2),value=TRUE))
expData2  <- as.matrix(getExpData(sN2))
res2      <- pruneKnn(expData2,no_cores=5)
cl2    <- graphCluster(res2,pvalue=0.01,use.leiden=TRUE, leiden.resolution=1.5)
#Noise estimation and related quantities:
x <- getFilteredCounts(sN2)
noise2 <- compTBNoise(res2,x,pvalue=0.01,gamma=1,no_cores=5)
qn2 <- quantKnn(res2, noise2, sN, pvalue = 0.01, minN = 5, no_cores = 5)

sN2 <- updateSC(sN2,res=res2,cl=cl2,noise=noise2,flo=.1)
sN2 <- comptsne(sN2)
sN2 <- compumap(sN2)

#UMAP plot
plotmap(sN2, um=TRUE,cex=.1)
#Plot of marker genes
fractDotPlot(sN2, genes, cluster=c(15,3,6,2,22,12,5,8,13,11,9,16,21,1,7,10,4,14,19,17,18,20), logscale=TRUE)
#Plot noise levels per cluster
clus2 <- c(15,3,6,2,22,12,5,13,8,11,9,16,21,1,7,10,4,14,19,17,18,20)
plotB(qn2$noise.av,sN2@cpart, cluster=3, set=clus2, pch = 20, cex = 0.2,col = sN2@fcol[clus2], ylab="biological noise",
      ylim=c(0,0.2))


#Perform batch correction and cell annotation in order to compare noise levels across similar cell populations
#Use harmony implemented within Seurat method
#Create a Seurat object
gen <- intersect(rownames(d), rownames(d2))
#Extract cell barcodes from the previous SCseq objects
cell1 <- colnames(sN@ndata)
cell2 <- colnames(sN2@ndata)
pbmcH <- CreateSeuratObject(counts = cbind(d[gen,cell1], d2[gen,cell1]), project = "PBMC")
pbmcH@meta.data$batch <- c(rep("Nuclei", length(cell1)), rep("WCell", length(cell2)))

#Seurat clustering
#Filtering
pbmcH <- subset(x = pbmcH, subset = nCount_RNA < 25000 & nCount_RNA > 1000)

#Initialize Seurat object
pbmcH <- Seurat::NormalizeData(pbmcH, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)
#Run harmony
pbmcH <- RunHarmony(pbmcH, "batch", plot_convergence = TRUE)
#Diagnostics plots
DimPlot(object = pbmcH, reduction = "harmony", pt.size = .1, group.by = "batch")
VlnPlot(object = pbmcH, features = "harmony_1", group.by = "batch", pt.size = .1)
#Dimension reduction, kNN network and clustering with harmony embeddings
pbmcH <- RunUMAP(pbmcH, reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
#UMAP plots with batch and clusters
DimPlot(pbmcH, reduction = "umap", group.by = "batch", pt.size = .1)
DimPlot(pbmcH, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, label=TRUE)

#Check contribution of batches into each cell type
table(pbmcH@meta.data$seurat_clusters)
df <- data.frame(Batch = pbmcH@meta.data$batch, Cluster = pbmcH@meta.data$seurat_clusters)
x <- paste(c("Nuclei", "WCell"), rep(0:18, each=2), sep="Clust")
df$BatchClust <- factor(paste(df$Batch, df$Cluster, sep="Clust"), levels=x)
head(df)
#Make a barplot
cols <- scales::hue_pal()(19)
barplot(table(df$BatchClust), col=rep(cols ,each=2), border=rep(c("black", "NA"), times=19),
        space=rep(c(0.25, 0.05), times=19))
#Cell type annotation based on marker genes
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmcH.markers <- FindAllMarkers(pbmcH, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
genesdf <- pbmcH.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
print(as.data.frame(genesdf))

#Make also expression UMAPs with gene expression profiles
#Chose one set at a time
FeaturePlot(pbmcH, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
FeaturePlot(pbmcH, features = c("IL7R", "CCR7", "CD14", "LYZ", "S100A4", "MS4A1", "CD8A", "FCGR3A", "MS4A7"))
FeaturePlot(pbmcH, features = c("GNLY", "NKG7", "FCER1A", "CST3", "PPBP", "PF4", "GP9"))

#Make a data frame in order to compare noise levels across cell populations
#Select cluster with unambiguous identity
clus <- c(4,2,0,6,9,8,1,3,11,10,5,7)
df <- data.frame(CellNames = rownames (pbmcH@meta.data[pbmcH@meta.data$seurat_clusters %in% clus,]))
df$Batch <- pbmcH@meta.data[df$CellNames,"batch"]
df$Batch[df$Batch == "WCell"] <- "Cell"
df$SeurClust <- pbmcH@meta.data[df$CellNames,"seurat_clusters"]
df$CellAnnot <- NA
df$CellAnnot[df$SeurClust == 0] <- "CD4 Memory"
df$CellAnnot[df$SeurClust %in% c(1,3,11)] <- "CD14 Monocytes"
df$CellAnnot[df$SeurClust == 2] <- "CD4 Naïve"
df$CellAnnot[df$SeurClust == 4] <- "CD8 Naïve"
df$CellAnnot[df$SeurClust %in% c(5,7)] <- "B cells"
df$CellAnnot[df$SeurClust == 6] <- "CD8 TEM1"
df$CellAnnot[df$SeurClust == 8] <- "NK"
df$CellAnnot[df$SeurClust == 9] <- "CD8 TEM2"
df$CellAnnot[df$SeurClust == 10] <- "CD16 Monocytes"
table(df$CellAnnot); sum(is.na(df$CellAnnot))


#Estimate average noise per cell across genes found in both datasets
k <- intersect(rownames(noise$epsilon), rownames(noise2$epsilon))
xe1 <- colMeans(noise$epsilon[k,], na.rm=TRUE)
xe2 <- colMeans(noise2$epsilon[k,], na.rm=TRUE)

df$Cell.noise <- NA
f <- df$CellNames[df$Batch == "Nuclei"]
x <- xe1[f]
df$Cell.noise[df$Batch == "Nuclei"] <- x

f2 <- df$CellNames[df$Batch == "Cell"]
x2 <- xe2[f2]
df$Cell.noise[df$Batch == "Cell"] <- x2
head(df)

#Violin plot
ggplot(df, aes(CellAnnot, Cell.noise, fill=Batch)) + geom_violin() +
  scale_fill_manual(values=c("red", "blue")) + ylab("Cellular Noise, caped to 0.4") + xlab("Cell types") +
  theme_classic() + ylim(0, 0.4)

#Scatter plot with error bars
#Take the information from the previous dataframe df and compute and stanrdard deviation per cell type and batch
x <- unique(df$CellAnnot)
dfNC <- data.frame(CellAnnot = x, MeanNoiseNuc = NA, sdNuc = NA, MeanNoiseCell = NA, sdCell = NA)
for (i in 1:length(x)){
  f <- df$Batch == "Nuclei" & df$CellAnnot %in% x[i]
  dfNC$MeanNoiseNuc[i] <- mean(df$Cell.noise[f], na.rm=TRUE)
  dfNC$sdNuc[i] <- sd(df$Cell.noise[f], na.rm=TRUE)
}
for (i in 1:length(x)){
  f <- df$Batch == "Cell" & df$CellAnnot %in% x[i]
  dfNC$MeanNoiseCell[i] <- mean(df$Cell.noise[f], na.rm=TRUE)
  dfNC$sdCell[i] <- sd(df$Cell.noise[f], na.rm=TRUE)
}

head(dfNC)
#Set some thresholds for the plot and colors
x2 <- c(dfNC[,2] + dfNC[,3], dfNC[,2] - dfNC[,3]) %>% range()
x3 <- c(dfNC[,4] + dfNC[,5], dfNC[,4] - dfNC[,5]) %>% range()
ax_lims <- range(c(x2, x3))
cols <- sN@fcol[1:9]
cols[5] <- "orange"
cols[7] <- "#255957"

plot(dfNC[,2], dfNC[,4], type="n", xlim=ax_lims, ylim=ax_lims, xlab="Cellular noise in Nuclei",
     ylab="Cellular noise in Cell")
arrows(dfNC[,2] - dfNC[,3], dfNC[,4], dfNC[,2] + dfNC[,3], dfNC[,4], length=0.05, angle=90, code=3, col=cols, lty=1, lwd=1.5)
arrows(dfNC[,2], dfNC[,4] - dfNC[,5], dfNC[,2], dfNC[,4] + dfNC[,5], length=0.05, angle=90, code=3, col=cols, lty=1, lwd=1.5)
points(dfNC[,2], dfNC[,4], pch=20, col=cols, cex=2)
abline(0,1, lty=2)
legend("topleft", dfNC[,1], col=cols, bty="n", pch=20)

#Focus in one cell populations and perform more detailed analyses:
#Select cells from CD8 Naive T cell compartment, found in cluster 4 of the pbmcH Seurat object
#TNa: T Naive cells
TNaNuc <- rownames (pbmcH@meta.data[pbmcH@meta.data$seurat_clusters == 4 & pbmcH@meta.data$batch == "Nuclei",])
TNaCell <- rownames (pbmcH@meta.data[pbmcH@meta.data$seurat_clusters == 4 & pbmcH@meta.data$batch == "WCell",])

#Perform differential expression analysis and select genes without significant expression change
k <- intersect(sN@genes,sN2@genes)
e1 <- sN@expdata[k , TNaNuc]
e2 <- sN2@expdata[k , TNaCell]
g <- rowSums(e1) + rowSums(e2) > 0
k <- k[g]
e1 <- e1[k,]
e2 <- e2[k,]
dgNC <- diffexpnb(cbind(e1,e2),A=colnames(e2),B=colnames(e1),method="per-condition", norm=TRUE)
#Extracting genes with differential expression
head(dgNC$res)
upNuc <- rownames(dgNC$res[dgNC$res$foldChange > 1.25 & dgNC$res$padj < 0.001,])
dnNuc <- rownames(dgNC$res[dgNC$res$foldChange < 0.8 & dgNC$res$padj < 0.001,])
#Genes without differential expression
geneF <- setdiff(k,c(upNuc,dnNuc))

#Ma plot for visualization:
mx <- dgNC$res$baseMean
names(mx) <- rownames(dgNC$res)
y <- dgNC$res$log2FoldChange
names(y) <- rownames(dgNC$res)
#Generate quantiles of the genes, based on their mean expression
quants <- quantile(mx, probs = seq(from=0, to=1, by=0.1), na.rm = TRUE)

plot(log2(mx),y,pch=20,cex=1,col="gray",xlab="log2 ((mean Nuclei + mean Cell)/2)",
     ylab="log2 (mean Nuclei) - log2 (mean Cell)")
points(log2(mx[geneF]), a9[geneF], pch=20, cex=1, col="red")
abline(v=log2(q[2:10]), col="black", lty=2, lwd=1.5)
legend("topright", c("Differentially Expressed", "No change in Expression"), pch=20, col=c("gray", "red"), bty="n")

#Take the no-differentialy expressed genes and compare noise values per each bin of genes, defined with the quants vector
#Aggregate mean expression values into bins:
x <- cut(mx[geneF], quants, right=FALSE, labels=FALSE) 
#Compute average noise per gene in the Naive T cells
#Note: Adding a pseudocount of 0.001
nm1 <- rowMeans(noise$epsilon[k , TNaNuc],na.rm=TRUE) + .001
nm2 <- rowMeans(noise2$epsilon[k , TNaCell],na.rm=TRUE) + .001

#Make a dataframe
df <- data.frame(MeanExpr = mx[geneF],
                 meanEpsNuc = nm1[geneF], meanEpsCell = nm2[geneF], BinExpr=x)
#df <- data.frame(MeanExprBoth = mx[a8], MeanExprNuc = m1[a8], meanExprCell = m2[a8], meanEpsNuc = nm1[a8], meanEpsCell = nm2[a8], BinExpr=a2)
head(df)
nL <- list()
j <- 1
for (i in 1:10){
  nL[[j]] <- df$meanEpsNuc [df$BinExpr == i]
  j <- j + 1
  nL[[j]] <- df$meanEpsCell [df$BinExpr == i]
  j <- j + 1
}
names(nL) <- paste( rep(1:10, each=2), rep(c("Nuclei", "Cell"), times=10), sep=":")
#Transform into log scale
nL2 <- lapply(nL, log2)
boxplot(nL2, col="white", border=rep(c("black", "red"), times=10), pch=20, cex=0.2, ylab="log2 (mean noise per gene)")

#Make a test for differential noisy genes between nuclear and cell transcripts
#Select CD8 Naive T cell, and common genes
n1 <- noise$epsilon[k , TNaNuc]
n2 <- noise2$epsilon[k , TNaCell]
#compute log2 fold change and p value
ln <- log2( rowMeans(n1,na.rm=TRUE)/rowMeans(n2,na.rm=TRUE) )
npv <-  apply(cbind(n1,n2), 1, function(x){
  y <- x[1:ncol(n1)]
  z <- x[(ncol(n1) + 1):length(x)]
  pval <- wilcox.test(z,y)$p.value
  return(pval)
  })
#Determine which are the differential noisy genes
f <- !is.na(ln) & abs(ln) > 1 & p.adjust(npv) < 1e-3
difNoiseG <- names(ln)[f]
#Have a look on some genes whose expression does not change, but noise changes
x <- ln[f]
x <- x[names(x) %in% geneF]
sort(x, decreasing=TRUE) %>% head(., 50)

#Make a MA plot for visualization
plot(log2(mx),ln,pch=20,col="grey", xlab="log2 ( (mean expression Nuclei + mean expression Cell) / 2 )",
     ylab="log2 (mean noise Nuclei / mean noise Cell)")
f <- !is.na(ln) & abs(ln) > 1 & p.adjust(npv) < 1e-3
points(log2(mx)[f],ln[f],pch=20,col="red")

#Selection of genes for single-molecule FISH experiemnts
#Check expression and noise patterns of gene candidates: PDCD4 and PPP1R2
plotexpmap(sN, "PDCD4", um=TRUE)
plotexpmap(sN, "PDCD4", um=TRUE, noise=TRUE)

plotexpmap(sN2, "PDCD4", um=TRUE)
plotexpmap(sN2, "PDCD4", um=TRUE, noise=TRUE)

plotexpmap(sN, "PPP1R2", um=TRUE)
plotexpmap(sN, "PPP1R2", um=TRUE, noise=TRUE)

plotexpmap(sN2, "PPP1R2", um=TRUE)
plotexpmap(sN2, "PPP1R2", um=TRUE, noise=TRUE)

