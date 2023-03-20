# This script contains the code to reproduce the analysis of aged and young hematopoietic stem cells from murine bone marrow with VarID2
# Contained in Figures 5, S5, 6a-b of Rosales-Alvarez et al., manuscript
# Data from Hérault, L., et al. Single-cell RNA-seq reveals a concomitant delay in differentiation and cell cycle of aged hematopoietic stem cells. BMC Biol. 2021 191 19, 1–20.
# The code in this script was performed with RaceID v0.2.6. For further versions, cluster numbers and parameter settings potentially change.

# Load required libraries
library(RaceID)
library(Matrix)
library(pheatmap)
library(ggplot2)
library(scales)

#Importing files
D <- list()
prefix <- c("/path_to_raw_data/")
n <- c("GSM4443875_young_A","GSM4443876_young_B","GSM4443877_old_A","GSM4443878_old_B")
for ( i in n ){
  d <- readMM(paste(prefix,i,"_matrix.mtx.gz",sep=""))
  f <- read.csv(paste(prefix,i,"_genes.tsv.gz",sep=""),sep="",header=FALSE)
  b <- read.csv(paste(prefix,i,"_barcodes.tsv.gz",sep=""),sep="",header=FALSE)
  
  j <- sub("GSM\\d+_","",i)
  rownames(d) <- make.unique(f[,2], sep = "_")
  colnames(d) <- paste(j,b[,1],sep="_")
  D[[j]] <- d
}
#Feature matrix with old cells
dO <- cbind(D[["old_A"]],D[["old_B"]])
#Feature matrix with old cells
dY <- cbind(D[["young_A"]],D[["young_B"]])
#Whole matrix
dOY <- cbind(dO,dY)
#Remove redundant objects
rm(D, dO, dY)


#Create a SCeq object and start with the standard analysis: quality filtering, knn inference and pruning, clustering and estimation of transition probabilities
sc <- SCseq(dOY)
sc <- filterdata(sc,mintotal=1000,CGenes=grep("^(mt|Gm\\d|Rp(l|s))",rownames(dOY),value=TRUE))
expData  <- getExpData(sc)
res      <- pruneKnn(expData,knn=25,no_cores=5)
cl    <- graphCluster(res,pvalue=0.01,min.size=5,use.leiden=FALSE)
probs <- transitionProbs(res,cl)

#Noise estimation and related quantities:
x <- getFilteredCounts(sc,minexpr=3,minnumber=5, gamma=0.5)
noise <- compTBNoise(res,x,pvalue=0.01,gamma=0.5,no_cores=5)
qn <- quantKnn(res, noise, sc, pvalue = 0.01, minN = 5, no_cores = 5)

#Add results to the SCseq object
sc <- updateSC(sc,res=res,cl=cl,noise=noise)

#Compute tSNE and UMAP representations
sc <- comptsne(sc,perplexity=200)
sc <- compumap(sc,n_neighbors=50,min_dist=1.5,spread=2)
plotmap(sc)
#Color by batches and age:
types <- sub("_[ACTGN]+\\-\\d+$","",colnames(sc@ndata))
plotsymbolsmap(sc,types, cex=0.3)

#Dotplot with marker genes
#Genes of interest
genes <- c("Slamf1","Hlf","Hoxa9","Mecom","Gata2","Kit","Ly6a","Cd48","Cd34","Cd74","Jun","Fos","Mki67","Flt3","Spi1","Vwf","Fli1")
#Order clusters in the desired way
clust_order <- c(7,15,1,6,3,4,5,9,13,10,14)
fractDotPlot(sc,genes=genes,cl=clust_order,zsc=TRUE,cap=0.5,flo=-0.5)

#Compute cellular noisex
qn <- quantKnn(res, noise, sc, pvalue = 0.01, minN = 5, map = "umap", iter = 10, no_cores = 5)

#Heatmap with the tSNE representation
plotQuantMap(qn,"noise.av",sc,box=FALSE,logsc=FALSE,ceil=.07)

#Comparison of cellular noise across clusters
#Clusters whose median would be highlighted
young_cl <- c(7,15)
#Boxplot with all clusters
plotQuantMap(qn,"noise.av",sc,box=TRUE,cluster=young_cl,logsc=TRUE, ylim=c(-6, -3))
#Only clusters corresponding to LT-HSCs
set <- c(7,15,1,6)
plotB(qn$noise.av, sc@cpart,cluster=young_cl,set=set,ylab="biological noise", col=sc@fcol[set], cex=0.2,ylim=c(0.01, 0.04))

#Wilcoxon test for the groups of young and aged LT-HSCs per batch:
a <- split(qn$noise.av, sc@cpart)
#Batch A: cluster 1 and 7
wilcox.test(a$"1", a$"7")
#Batch B: cluster 6 and 15
wilcox.test(a$"6", a$"15")

#Differential noise analysis between old and young LT-HSCs:
#Clusters with old LT-HSCs:
setO <- c(1,6)
#Clusters with young LT-HSCs:
setY <- c(7,15)
#Perform the test:
ngenes <- diffNoisyGenesTB(noise,cl,setO,setY,no_cores=5)
plotDiffNoise(ngenes, lthr=log2(1.25), show_names=TRUE)


#Check noise levels for some genes of interest
g <- c("Dlk1", "Terf1", "Gfi1b", "Cyp26b1", "Cdkn2c")
#Create a data frame with noise log-transformed noise signal
for ( i in 1:length(g) ){
  d <- noise$epsilon
  
  x1 <- d[g[i],sc@cpart %in% setY[1]]
  x2 <- d[g[i],sc@cpart %in% setY[2]]
  y1 <- d[g[i],sc@cpart %in% setO[1]]
  y2 <- d[g[i],sc@cpart %in% setO[2]]
  
  if (i == 1)
    df <- data.frame( y = log( c( x1, x2, y1, y2) + .05 ) , x = g[i], Group = c( rep( "young", length( c(x1,x2) ) ),  rep( "old", length( c(y1,y2) ) ) ) )
  else
  {df2 <- data.frame( y = log( c( x1, x2, y1, y2) + .05 ) , x = g[i], Group = c( rep( "young", length( c(x1,x2) ) ),  rep( "old", length( c(y1,y2) ) ) ) ) 
  df <- rbind (df, df2)
  }
}

library(ggplot2)
p <- ggplot(df, aes(x, y, fill=Group)) + geom_violin(scale="width") + scale_fill_manual(values=c("red", "blue")) +
  ylab("noise") + xlab("Gene") + ggtitle("Noise levels") + ylim(-3,0) + theme_classic()


#Check expression of Dlk1 in the tSNE plot (related to Figure 6A in the manuscript)
plotexpmap(sc, "Dlk1", log=TRUE,noise=FALSE,cex=.75)
#Check also noise of Dlk1 in the tSNE plot
plotexpmap(sc, "Dlk1", log=TRUE,noise=TRUE,cex=.75)


#Differential expression analysis of cells expressing Dlk1 or not expressing it (related to Figure 6B in the manuscript)
#Comparison within the LT-HSCs from old mice
#Old samples, batch A
oldA <- 1
x <- getFilteredCounts(sc,minexpr=3,minnumber=5)
x1 <- x[,cl$partition %in% oldA]
A <- colnames(x1)[x1["Dlk1",] == 0]
B <- colnames(x1)[x1["Dlk1",] > 0]
z1 <- diffexpnb(x1,A,B,norm=TRUE)
plotdiffgenesnb(z1)

#Old samples, batch B
oldB <- 6
x2 <- x[,cl$partition %in% oldB]
A <- colnames(x2)[x2["Dlk1",] == 0]
B <- colnames(x2)[x2["Dlk1",]> 0]
z2 <- diffexpnb(x2,A,B,norm=TRUE)
plotdiffgenesnb(z2)