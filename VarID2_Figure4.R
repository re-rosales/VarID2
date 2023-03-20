# This script contains the code to reproduce the analysis of hematopoietic stem cells with VarID2
# Contained in Figures 4 and S4 of Rosales-Alvarez et al., manuscript
# Data from Dahlin,  et al. A single-cell hematopoietic landscape resolves 8 lineage trajectories and defects in Kit mutant mice. Blood. 2018;131:e1–11.
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
library(FateID)
library(quadprog)
library(pheatmap)
library(ReactomePA)
library(clusterProfiler)

#Import files
n <- c("27_SIGAB1","28_SIGAC1","29_SIGAD1","30_SIGAF1","31_SIGAG1","32_SIGAH1","33_SIGAG8","34_SIGAH8")
j <- 0
for ( i in n ){
  x <- Matrix(as.matrix(read.csv(paste("path_to_data/",i,"_counts.txt",sep=""),header=TRUE,sep="\t")),sparse=TRUE)
  colnames(x) <- paste(i,colnames(x),sep=".")
  if ( j == 0 ){
    j <- 1
    D <- x
  }else{
    D <- cbind(D,x)
  }
}

#Arrange rownames
library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- rownames(D)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)
map <- merge(data.frame(gene=genes),G_list,by.x="gene",by.y="ensembl_gene_id",all.x=TRUE)
map <- t(apply(map,1,function(x){ if ( is.na(x[2]) ){ x[2] <- x[1]}; return(x) }))
rownames(map) <- map[,1]
map <- map[rownames(D),]
rownames(D) <- make.unique(map[,2], sep = "_")

#Analysis of wild type samples, contained in the first 6 files:
#This a large dataset and the following steps (control quality, knn inference, pruning, clustering, UMAP and TSNE representations, and noise inference) can take a long time to run (a couple of hours).
#We suggested to run it on batch mode through RScript
sw <- SCseq(D[,grep("33_SIGAG8|34_SIGAH8",colnames(D),invert=TRUE)])
sw <- filterdata(sw,mintotal=2000,CGenes=grep("^(mt|Gm\\d|Rp(l|s))",rownames(D),value=TRUE))
expData <- getExpData(sw)
resw     <- pruneKnn(expData,no_cores=5,seed=12345,knn=50,pcaComp=30)
clw      <- graphCluster(resw,pvalue=0.01,use.leiden=TRUE,leiden.resolution=1.5)
probsw   <- transitionProbs(resw,clw)
sw <- updateSC(sw,res=resw,cl=clw)
sw <- comptsne(sw,perplexity=200)
sw <- compumap(sw,min_dist=1.5,spread=2)

plotmap(sw, um=TRUE, cex=0.3)
plotTrProbs(sw,probsw,tp=.5,prthr=0.01,cthr=0,um=TRUE, cex=0.3)

#Dot plot for expression of some relevant marker genes
genes <- c("Slamf1","Hlf","Mecom","Hoxa9","Kit","Ly6a","Spi1","Cd34","Cd48","Jun","Fos","Cd74",
           "Mki67","Flt3","Mcpt8","Vwf","Fli1","Gata2","Gata1")
#Order clusters
clusw <- c(10,1,11,15,18,7,4,12,14,9,2,22,17,5,19,6,20,3,13,21,8,16,23)
fractDotPlot(sw,genes=genes,cl=clusw,zsc=TRUE, cap=0.75,flo=-0.75)

#Compute noise
x <- getFilteredCounts(sw,minexpr=5,minnumber=5)
noisew <- compTBNoise(resw,x,pvalue=.01,no_cores=5,gamma=0.5,x0=0,lower=0,upper=100)
sw <- updateSC(sw,res=resw,cl=clw,noise=noisew)
rm(x)
qw <- quantKnn(resw, noisew, sw, pvalue = 0.01, minN = 5, no_cores = 5)

#Exploring some output from noise analysis
#Noise levels across cell clusters
StemW <- 10
plotB(log2(qw$noise.av),sw@cpart, cluster=StemW, set=clusw, cex = 0.2, col = sw@fcol[clusw],
      ylab="biological noise", ylim=c(-5, -2))
#UMAP representations with UMI counts per cell
plotQuantMap(qw,"umi",sw,um=TRUE,logsc=TRUE, cex=0.3)
#UMAP representations with local cell-cell correlations
plotQuantMap(qw,"local.corr",sw,um=TRUE,logsc=TRUE, cex=0.3)
#UMAP representations with cellular noise levels
plotQuantMap(qw,"noise.av",sw,um=TRUE, logsc=TRUE, cex=0.3)

#Check correlations between cellular noise, UMI counts per cell and local cell-cell correlations
plotQQ(qw,"local.corr","noise.av",sw,cluster=StemW,log="xy")
plotQQ(qw,"umi","noise.av",sw,cluster=StemW,log="yx")
plotQQ(qw,"umi","local.corr",sw,cluster=StemW,log="yx")


#Analysis of the mutant samples
sk <- SCseq(D[,grep("33_SIGAG8|34_SIGAH8",colnames(D))])
sk <- filterdata(sk,mintotal=2000,CGenes=grep("^(mt|Gm\\d|Rp(l|s))",rownames(D),value=TRUE))
expData  <- getExpData(sk)
resk   <- pruneKnn(expData,no_cores=5,seed=12345,knn=50,pcaComp=30)
clk    <- graphCluster(resk,pvalue=0.01,use.leiden=TRUE,leiden.resolution=1.5)
probsk <- transitionProbs(resk,clk)
sk <- updateSC(sk,res=resk,cl=clk)
sk <- comptsne(sk,perplexity=200)
sk <- compumap(sk,min_dist=1,spread=1.5)

plotmap(sk,um=TRUE, cex=0.3)
plotTrProbs(sk,probsw,tp=.5,prthr=0.01,cthr=0,um=TRUE, cex=0.3)
#Dot plot for expression of some relevant marker genes
genes <- c("Slamf1","Hlf","Mecom","Hoxa9","Kit","Ly6a","Spi1","Cd34","Cd48","Jun","Fos","Cd74",
           "Mki67","Flt3","Mcpt8","Vwf","Fli1","Gata2","Gata1")
clusk <- c(17,1,7,15,18,20,10,4,16,13,5,21,6,8,19,11,22,9,12,3,2,14)
fractDotPlot(sk,genes=genes, cl=clusk, zsc=TRUE, cap=0.75,flo=-0.75)

#Noise analysis:
x <- getFilteredCounts(sk,minexpr=5,minnumber=5)
noisek <- compTBNoise(resk,x,pvalue=.01,no_cores=5,gamma=0.5,x0=0,lower=0,upper=100)
sk <- updateSC(sk,res=resk,cl=clk,noise=noisek)
rm(x)
qk <- quantKnn(resk, noisek, sk, pvalue = 0.01, minN = 5, no_cores = 5)

#Noise levels across cell clusters
StemK <- 17
clusk <- c(17,1,7,15,18,20,10,4,16,13,5,21,6,8,19,11,22,9,12,3,2,14)
plotB(log2(qk$noise.av),sk@cpart,cluster=StemK,ylab="biological noise", set=clusk, col=sk@fcol[clusk],
      ylim=c(-4.8,-3.2), cex=0.2)

#Other related quantities
plotQuantMap(qk,"umi",sk,um=TRUE,logsc=TRUE)
plotQuantMap(qk,"local.corr",sk,um=TRUE,logsc=TRUE)
plotQuantMap(qk,"noise.av",sk,um=TRUE,logsc=TRUE)
#Linear scale:
plotQuantMap(qk,"noise.av",sk,um=TRUE, ceil=0.1)


#Quadratic programming to match similar cell types
QP <- function(k,m,norm=TRUE){
  Dmat <- t(m) %*% m
  dvec <- t(k) %*% m
  if ( norm ){
    Amat <- cbind(rep(1,ncol(m)), diag(ncol(m)))
    bvec <- c(1,rep(0,ncol(m)))
    qp <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  }else{
    Amat <- diag(ncol(m))
    bvec <- rep(0,ncol(m))
    qp <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0, factorized=FALSE)
  }
  
  return( list(w=qp$solution, fit=m %*% qp$solution, residual= sum((m %*% qp$solution - k)**2), qp=qp))
}

k <- intersect( sw@genes, sk@genes )

for ( i in 1:max(sk@cpart) ){
  if ( i  == 1 ) mk <- rowMeans(as.matrix(sk@ndata[k,sk@cpart == i])) else mk <- cbind(mk, rowMeans( as.matrix(sk@ndata[k,sk@cpart == i]))) 
}
colnames(mk) <- paste("K",1:ncol(mk),sep=".")

for ( i in 1:max(sw@cpart) ){
  if ( i  == 1 ) mw <- rowMeans( as.matrix(sw@ndata[k,sw@cpart == i]) ) else mw <- cbind(mw, rowMeans( as.matrix(sw@ndata[k,sw@cpart == i]) )) 
}
colnames(mw) <- paste("W",1:ncol(mw),sep=".")

f <- apply(mk,1,var) > 0 & apply(mw,1,var) > 0
m <- cor(as.matrix(mw[f,]),as.matrix(mk[f,]),use="pairwise.complete.obs")
m <- cor(as.matrix(mw[f,]),as.matrix(mk[f,]),use="pairwise.complete.obs",method="spearman")

w <- matrix(rep(NA,ncol(mw)*ncol(mk)),ncol=ncol(mw) )
for ( i in 1:ncol(mk) ){
  w[i,] <- QP(mk[f,i],mw[f,])$qp$solution   
}

colnames(w) <- colnames(mw)
rownames(w) <- colnames(mk)
#w contains the similarity weights of clusters in the W41/W41 dataset (in rows) to WT clusters (in columns) inferred by quadratic programming
#Plot for better visualization
pheatmap(w,cluster_cols=FALSE,cluster_rows=FALSE)


#Comparison between LT-HSCs in KO versus WT datasets, clusters 17 and 10, respectively
StemK <- 17
StemW <- 10

#Differential gene expression analysis
#Take expression data
k <- intersect( sw@genes, sk@genes )
ek <- as.matrix(sk@expdata[k , colnames(sk@ndata)[clk$partition %in% StemK]])
ew <- as.matrix(sw@expdata[k , colnames(sw@ndata)[clw$partition %in% StemW]])
#Filter for expressed gnees in these subsets
g <- rowSums(ek) + rowSums(ew) > 0
k <- k[g]
ek <- ek[k,]
ew <- ew[k,]

dg <- diffexpnb(cbind(ek,ew),A=colnames(ew),B=colnames(ek),norm=TRUE)
plotdiffgenesnb(dg,lthr=log2(1.25),show_names=TRUE,pthr=0.001, Aname="WT", Bname="KitMut")

#Differential noisy genes
nk <- noisek$epsilon[k , clk$partition %in% StemK]
nw <- noisew$epsilon[k , clw$partition %in% StemW]

#Log2 Fold change:
ln <- log2( rowMeans(nk,na.rm=TRUE)/rowMeans(nw,na.rm=TRUE) )
#p values:
npv <-  apply(cbind(nk,nw),1,function(x){ y <- x[1:ncol(nk)]; z <- x[(ncol(nk) + 1):length(x)]; wilcox.test(z,y)$p.value } )

#Take mean expression and mean noise for MA plots
mk <- rowMeans(ek)
mw <- rowMeans(ew)
ak <- rowMeans(noisek$epsilon[k , clk$partition %in% StemK],na.rm=TRUE) + .001
aw <- rowMeans(noisew$epsilon[k , clw$partition %in% StemW],na.rm=TRUE) + .001
plot((mk+mw)/2,ak/aw,pch=20,col="grey",log="xy")
f <- !is.na(ln) & abs(ln) > 1 & p.adjust(npv) < 1e-3
points((mk+mw)[f]/2,(ak/aw)[f],pch=20,col="red")
text((mk+mw)[f]/2,(ak/aw)[f],k[f],cex=.5)



plot((mk+mw)/2,ak/aw,pch=20,col="grey",log="xy")
f <- !is.na(ln) & abs(ln) > 1 & p.adjust(npv) < 1e-3
points((mk+mw)[f]/2,(ak/aw)[f],pch=20,col="red")
text((mk+mw)[f]/2,(ak/aw)[f],k[f],cex=.5)
#png("D10KOvsWT04_DiffVar02.pdf")
plot((mk+mw)/2,ak/aw,pch=20,col="grey",log="xy")
#f <- !is.na(ln) & abs(ln) > 1 & p.adjust(npv) < 1e-3 #repeated
points((mk+mw)[f]/2,(ak/aw)[f],pch=20,col="red")


#Check which genes increase expression:
#In mutant data
dgK <- rownames(dg$res)[dg$res$foldChange > 1.25 & dg$res$padj < 0.001]
length(dgK)
head(dgK, 50)
#In wild type data
dgW <- rownames(dg$res)[dg$res$foldChange < 0.8 & dg$res$padj < 0.001]
length(dgW); head(dgW, 300)

#Check which genes increase noise:
#In mutant data
dNK <- names(ln) [!is.na(ln) & ln > 1 & p.adjust(npv) < 1e-3]
length(dNK)
head(dNK, 100)
#In wild type data
dNW <- names(ln) [!is.na(ln) & ln < -1 & p.adjust(npv) < 1e-3]
length(dNW)
head(dNW, 100)

#Enrichment analysis of over-expressed genes in LT-HSC (dgK vector)
de <- dgK
de <- bitr(de, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
uni <- intersect(sw@genes, sk@genes)
uni <- bitr(uni, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

x <- enrichPathway(gene=de[,2], organism="mouse", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable=TRUE, universe=uni[,2])
dim(x)
#No enrichment

#Enrichment analysis of noisy genes in old qHSC (gVk2 vector)
de <- dNK
de <- bitr(de, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
uni <- intersect(sw@genes, sk@genes)
uni<-bitr(uni, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

x <- enrichPathway(gene=de[,2], organism="mouse", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable=TRUE, universe=uni[,2])
dim(x)
barplot(x, showCategory=11)

#Expression and noise profiles of some genes of interest
#Compute data frames and violin plots
g <- c("Pf4", "Mcm2", "Mcm3", "Mcm5", "Mcm7","Mpl","Ccl9","Orc6", "Cdt1")
df <- data.frame()

shex <- TRUE #Indicator to take gene expression values and make a data frame
for ( i in 1:length(g) ){
  if ( shex ){
    dk <- sk@ndata
    dw <- sw@ndata
    la <- "expression"
  }else{
    dk <- noisek$epsilon
    dw <- noisew$epsilon
    la <- "noise"
  }
  x1 <- dk[g[i],sk@cpart == StemK]
  x2 <- dk[g[i],sk@cpart != StemK]
  y1 <- dw[g[i],sw@cpart == StemW]
  y2 <- dw[g[i],sw@cpart != StemW]
  
  if (i == 1)
    df <- data.frame( y = log( c( x1, x2, y1, y2) + .05 ) , x = g[i], Sample =  c( rep("W41 HSC",length(x1)), rep("W41 MPP",length(x2)), rep("WT HSC",length(y1)), rep("WT MPP",length(y2)) ) )
  else {
    df2 <- data.frame( y = log( c( x1, x2, y1, y2) + .05 ) , x = g[i], Sample =  c( rep("W41 HSC",length(x1)), rep("W41 MPP",length(x2)), rep("WT HSC",length(y1)), rep("WT MPP",length(y2)) ) )
    df <- rbind (df, df2)
  }
}

gg <- ggplot(df, aes(x, y, fill=Sample)) + geom_violin() +
  scale_fill_manual(values=c("red", "#ff7f7f","blue","#7f7fff")) +
  ylab(la) + xlab("Gene") + ylim(min(df$y),-2.95) + theme_classic()
gg

#Repeat the same loop, this time for noise data:
df <- data.frame()
shex <- FALSE #Indicator to take noise values this time
for ( i in 1:length(g) ){
  if ( shex ){
    dk <- sk@ndata
    dw <- sw@ndata
    la <- "expression"
  }else{
    dk <- noisek$epsilon
    dw <- noisew$epsilon
    la <- "noise"
  }
  x1 <- dk[g[i],sk@cpart == StemK]
  x2 <- dk[g[i],sk@cpart != StemK]
  y1 <- dw[g[i],sw@cpart == StemW]
  y2 <- dw[g[i],sw@cpart != StemW]
  
  if (i == 1)
    df <- data.frame( y = log( c( x1, x2, y1, y2) + .05 ) , x = g[i], Sample =  c( rep("W41 HSC",length(x1)), rep("W41 MPP",length(x2)), rep("WT HSC",length(y1)), rep("WT MPP",length(y2)) ) )
  else {
    df2 <- data.frame( y = log( c( x1, x2, y1, y2) + .05 ) , x = g[i], Sample =  c( rep("W41 HSC",length(x1)), rep("W41 MPP",length(x2)), rep("WT HSC",length(y1)), rep("WT MPP",length(y2)) ) )
    df <- rbind (df, df2)
  }
}

gg2 <- ggplot(df, aes(x, y, fill=Sample)) + geom_violin(scale="width") +
  scale_fill_manual(values=c("red", "#ff7f7f","blue","#7f7fff")) +
  ylab(la) + xlab("Gene") + theme_classic()
gg2

