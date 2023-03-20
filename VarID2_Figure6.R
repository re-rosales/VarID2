# This script contains the code to reproduce the analysis of Dlk1+ and Dlk1- long-term hematopoietic stem cells (LT-HSCs)
# Contained Figures 6c-e and S6a-b of Rosales-Alvarez et al., manuscript
# Data generated with mCELseq2 protocol. Data available on GEO database, accession number GSE185637
# The code in this script was performed with RaceID v0.2.6. For further versions, cluster numbers and parameter settings potentially change.

# Setup python environment in order to use leiden clustering
# It requires to install reticulate R package and the python "leidenalg" and "igraph" modules in the indicated python environment
reticulate::use_python("./anaconda3/bin/python", required = TRUE) #if required, change the path of the corresponding environment environment
# Corroborate "leidenalg" and "igraph" are available
reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()


library(Matrix)
library(RaceID)
library(readr)
library(data.table)
library(dplyr)
library(tools)

#Import data
prefix <- "/path_to_CEL-Seq2_mapped_files/"
files <- list.files(path=prefix)
files <- grep("coutt.csv", files, value=TRUE)
#Read data
data_list <- lapply(paste(prefix, files, sep="/"), function(x) {fread(x, header= T)} )
names(data_list) <- sub("\\.coutt.csv", "",files)
#Make a reference  data frame with all the common gene names
geneids <- lapply(data_list, function(x) x$GENEID)
geneid_ref <- unlist(geneids) %>% unique()
geneid_ref <- data.frame("GENEID"=geneid_ref)

#Assign column names and add all the gene names found
for (d in names(data_list)) { 
  colnames(data_list[[d]]) <- c("GENEID", paste(d, "_",1:192,sep="" )); #192 is the number of cells in each plate, corresponding to ncol without GENEID
  setkey(data_list[[d]], GENEID); 
  data_list[[d]] <- merge( geneid_ref, data_list[[d]], by="GENEID",all.x = T )  
}

# Combine list of data.tables and removing the GENEID column from data.tables
data_list_cbind <- dplyr::bind_cols( lapply(data_list, function(x) { x[ order(match(x$GENEID, geneid_ref$GENEID)), ][-1] }) )
#Change NA values into 0
data_list_cbind[is.na(data_list_cbind)] <- 0

#Transform into a matrix and assign row names
data <- as.matrix(data_list_cbind)
rownames(data) <- geneid_ref$GENEID
#Remove no-necessary objects
rm(data_list, data_list_cbind, geneids, geneid_ref)

#Identify contaminating red blood cells with high expression of Hba-a1
x <- data["Hba-a1",]
cont <- names(x[x >10])
f <- setdiff(colnames(data), cont)
data <- data[,f]

#Create a SCeq object, perform quality filtering, knn inference and pruning, clustering and dimension reduced representations
sC <- SCseq(data)
sC <- filterdata(sC,mintotal=2000,CGenes=grep("^(mt|Gm\\d|Rp(l|s))",rownames(data),value=TRUE), minexpr = 3,
                  minnumber = 3)
eD <- getExpData(sC)
resc     <- pruneKnn(eD,no_cores=5,seed=12345)
clc      <- graphCluster(resc2pvalue=0.01,use.leiden=TRUE,leiden.resolution=1)
probsc   <- transitionProbs(resc,clc)
sC <- updateSC(sC,res=resc,cl=clc)
sC <- comptsne(sC)
sC <- compumap(sC, min_dist=2,spread=3)

#UMAP plot
plotmap(sC, um=TRUE)
#UMAP with Dlk1 positive and negative cells
f <- sub("__.+","",colnames(eD))
plotsymbolsmap(sC, f, um=TRUE)
#UMAP with DLk1 expression
plotexpmap(sC2, "Dlk1", um=TRUE, logsc=TRUE, cex=1)

####Symbols map, Dlk1 expression, dot plot and differential expression analysis

#Dot plot with expression of relevant genes
genes <- c("Jun","Fos","Klf6","Nr4a1","Mecom","Hlf","Ly6a","Procr", "Hoxa9", "Kit","Dlk1", "Meg3","Sult1a1", 
        "Mki67", "Ccnb1", "Pcna","Cd48","Cd34")
clusterV <- c(1,4,2,3,5)
fractDotPlot(sC, genes, cluster=clusterV,zsc=TRUE,cap=1)

#Perform differential expression analysis between Dlk1 positive and negative cells
#f <- sub("__.+$","", colnames(eD))
A <- colnames(eD)[grepl("neg_old", colnames(eD))]
B <- colnames(eD)[grepl("pos_old", colnames(eD))]
diffG <- diffexpnb(getfdata(sC,n=c(A,B)), A=A, B=B )

#Plot with differential expressed genes
#Threshold values: log2FC > log2(1.25), padj < 0.05
plotdiffgenesnb(diffG,pthr=.05,lthr=log2(1.25),show_names=TRUE,padj=TRUE)
