# This script contains the code to reproduce the analysis of the multiomics data with VarID2
# Contained in Figures 3 and 3S of Rosales-Alvarez et al., manuscript
# Data from 10X genomics, samples from peripheral blood mononuclear cells (PBMC) from a healthy donor - granulocytes removed through cell sorting (10k)
# The code in this script was performed with RaceID v0.2.6. For further versions, cluster numbers and parameter settings potentially change.

# Setup python environment in order to use leiden clustering
# It requires to install reticulate R package and the python "leidenalg" and "igraph" modules in the indicated python environment
reticulate::use_python("./anaconda3/bin/python", required = TRUE)
# Corroborate "leidenalg" and "igraph" are available, as well as python
reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()

# Load required libraries
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)
library(RaceID)
library(Matrix)
library(pheatmap)
library(scales)
library(magrittr)
library(ReactomePA)
library(clusterProfiler)
library(RcisTarget)

#Import multiomics data
#Download input data from: 10x Genomics,
#https: //support.10xgenomics.com/single-cell-multiome-atac-gex/ datasets/1.0.0/pbmc_granulocyte_sorted_10k
inputdata.10x <- Read10X_h5("/path_to_data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")

#Follow Signac vignette
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#Add the ATAC-seq data
#Use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/path_to_data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

#Basic QC based on the number of detected molecules for each modality as well as mitochondrial percentage.
VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

#Filtering
dim(expData)
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

#pre-processing and dimensional reduction on both assays independently, using standard approaches for RNA and
#ATAC-seq data
# RNA analysis
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. This graph is further
#used for UMAP visualization and clustering
pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# perform sub-clustering on cluster 6 to find additional structure
pbmc <- FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(pbmc) <- "sub.cluster"

# add annotations
pbmc <- RenameIdents(pbmc, '19' = 'pDC','20' = 'HSPC','15' = 'cDC')
pbmc <- RenameIdents(pbmc, '0' = 'CD14 Mono', '9' ='CD14 Mono', '5' = 'CD16 Mono') #0 and 9 are joined
pbmc <- RenameIdents(pbmc, '10' = 'Naive B', '11' = 'Intermediate B', '17' = 'Memory B', '21' = 'Plasma')
pbmc <- RenameIdents(pbmc, '7' = 'NK')
pbmc <- RenameIdents(pbmc, '4' = 'CD4 TCM', '13'= "CD4 TEM", '3' = "CD4 TCM", '16' ="Treg", '1' ="CD4 Naive", '14' = "CD4 Naive")
pbmc <- RenameIdents(pbmc, '2' = 'CD8 Naive', '8'= "CD8 Naive", '12' = 'CD8 TEM_1', '6_0' = 'CD8 TEM_2', '6_1' ='CD8 TEM_2', '6_4' ='CD8 TEM_2')
pbmc <- RenameIdents(pbmc, '18' = 'MAIT')
pbmc <- RenameIdents(pbmc, '6_2' ='gdT', '6_3' = 'gdT')
pbmc$celltype <- Idents(pbmc)

#visualize clustering based on gene expression, ATAC-seq, or WNN analysis
p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

#Additional plot:
DimPlot(pbmc, reduction = "wnn.umap", group.by = "seurat_clusters", pt.size = .1, label=TRUE)

#Compute markers genes for each cluster
DefaultAssay(pbmc) <- "RNA"
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(pbmc.markers)
genes <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
print(as.data.frame(genes))




#Compute gene activities
DefaultAssay(pbmc) <- "ATAC"
gene.ac <- GeneActivity(pbmc, max.width = NULL)

#Normalize by kilobase gene length and then by mean counts per cell
#Generate a genomic ranges object with the length of genes and 2kb upstream of the transcriptional start site
#Start by extracting genomic locations with biomaRt
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
transcripts <- getBM(attributes=c('chromosome_name', 'start_position','end_position','strand', 'ensembl_gene_id',
                                  "external_gene_name", "gene_biotype"), filters ="biotype",
                     values  =c("protein_coding"), mart =mart)
head(transcripts)
#Build a GRanges with regioneR's toGRanges function
transcripts.gr <- toGRanges(transcripts)
#Filter out non-standard chromosomes
transcripts.gr <- transcripts.gr[seqnames(transcripts.gr) %in% c(1:22, "X", "Y")]
#Add strand sign, based on the metadata slot of the same object
strand(transcripts.gr) [transcripts.gr$strand == 1] <- "+"
strand(transcripts.gr) [transcripts.gr$strand == -1] <- "-"
#Once having strand sign, remove the strand column of the metadata
transcripts.gr <- transcripts.gr[,c("ensembl_gene_id", "external_gene_name", "gene_biotype")]
#Expand 2000bp upstream to the genomic ranges
transcripts.gr <- promoters(transcripts.gr, upstream=2000, downstream = width(transcripts.gr))
#Set to USCS style and genome = hg38
seqlevelsStyle(transcripts.gr) <- 'UCSC'
genome(transcripts.gr) <- "hg38"
transcripts.gr

#Extract the width of each gene
x <- width(transcripts.gr)
hist(a, breaks=100)
names(x) <- transcripts.gr$external_gene_name
head(x)
x <- x/1000 #scale to 1000 bp

#Select overlapped genes
k <- intersect(rownames(gene.ac), transcripts.gr$external_gene_name)
length(k)
x <- x[k]
sum(duplicated(names(x)))

#Normalize gene activities by kilobase length
gene.ac2 <- gene.ac[k,] / x[k]
dim(gene.ac2)

#Normalize by mean counts per cell based on the activity matrix
mx <- colMeans(gene.ac[k,])
gene.ac2 <- gene.ac2 / mx
#Add a pseudo count for 0 values
gene.ac2[gene.ac2 == 0] <- 1e-05


#Compute noise with VarID2, similar as in Fig. 2:
d <- inputdata.10x$`Gene Expression`
#Create a SCseq object, perform quality filtering, knn inference and pruning, clustering and noise estimation
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
plotmap(sN, um=TRUE)

#Make differential expression analysis of the major cell types: T cells, Monocytes and B cells
difG <- list()
difG$"T cells" <- clustdiffgenes(sN, cl = c(16,2,17,3,11,18,1,14,20,6), bgr = NULL, pvalue = 0.01)
difG$"Mono" <- clustdiffgenes(sN, cl = c(4,7,5,9,15,8), bgr = NULL, pvalue = 0.01)
difG$"B cells" <- clustdiffgenes(sN, cl = c(10,13,22), bgr = NULL, pvalue = 0.01)

#For visualization of gene activity scores into RaceID/VarID functions, make a copy of the SCseq object
#Add the gene activities into the noise slot
sN2 <- sN
sN2@noise <- gene.ac2

#Compute correlations between gene expression and gene activities, and noise and gene activities
#Normalized gene expression:
nData <- getfdata(sN)
#Noise:
eps <- noise$epsilon
eps[is.na(eps)] <- 0
#Set intersect of genes
k <- intersect(rownames(gene.ac2), rownames(nData))
k <- intersect(k, rownames(eps))

#Correlations expression and gene activity
ag <- cbind(nData[k,], gene.ac2[k,])

agC <- apply(ag,1,function(x){
  x2 <- x[1:ncol(nData)]
  x3 <- x[(ncol(nData)+1):length(x)]
  x4 <- cor(x2,x3,use="pairwise.complete.obs", method="pearson")
  return(x4)
  })

#Quick check of scores
length(agC)
sum(agC > 0.05)
sum(agC < -0.05)

#Correlations noise and gene activity
ag <- cbind(eps[k,], gene.ac2[k,])

ag2C <- apply(ag,1,function(x){
  x2 <- x[1:ncol(eps)]
  x3 <- x[(ncol(eps) +1):length(x)]
  x4 <- cor(x2,x3,use="pairwise.complete.obs", method="pearson")
  return(x4)
})

length(ag2C)
sum(ag2C > 0.05)
sum(ag2C < -0.05)

#Check intersects, defining class A, C and B genes
Agen <- intersect(names(agC)[agC > 0.05] , names(ag2C)[ag2C > 0.05])
Cgen <- intersect(names(agC)[agC > 0.05] , names(ag2C)[ag2C < 0.05 & ag2C > -0.05])
Bgen <- intersect(names(agC)[agC > 0.05] , names(ag2C)[ag2C <-0.05])
length(Agen)
head(Agen, 50)
length(Bgen)
head(Bgen, 50)
length(Cgen)
head(Cgen, 50)

#Check distribution scores for class A and B genes
#Expression and gene activity correlations:
h <- hist(agC[c(Agen, Bgen)],breaks=40,plot=FALSE)
xa <- range(h$breaks)
h1 <- hist(agC[Agen],breaks=h$breaks,plot=FALSE)
h2 <- hist(agC[Bgen],breaks=h$breaks,plot=FALSE)
#Noise and gene activity correlations:
hn <- hist(ag2C[c(Agen, Bgen)],breaks=40,plot=FALSE)
xa2 <- range(hn$breaks)
h3 <- hist(ag2C[Agen],breaks=hn$breaks,plot=FALSE)
h4 <- hist(ag2C[Bgen],breaks=hn$breaks,plot=FALSE)

colsa <- c("#008744", "#9f3268")
plot(h1,border="white", col=scales::alpha(colsa[1], .4),freq=FALSE, xlim=xa,
     ylim=c(0,max(c(h1$density, h2$density))), xlab="Pearson Correlation scores",
     main="Correlation: Expression - Gene Activity")
plot(h2,border="white",col=scales::alpha(colsa[2], .4),freq=FALSE, add=TRUE)
legend("topright", legend=c("Class A", "Class B"), pch=15, col=scales::alpha(colsa, 0.4), bty="n")

plot(h3,border="white",col=scales::alpha(colsa[1], .4),freq=FALSE, xlim=xa2,
     ylim=c(0,max(c(h3$density, h4$density))), xlab="Pearson Correlation scores",
     main="Correlation: Noise - Gene Activity")
plot(h4,border="white",col=scales::alpha(colsa[2], .4),freq=FALSE, add=TRUE)
legend("topright", legend=c("Class A", "Class B"), pch=15, col=scales::alpha(colsa, 0.4), bty="n")


#Heat plots for checking patterns of expression, noise and gene activity
#Make a subset of class A genes
#Take results from differential expression analysis for T cells, B cells and Monocytes, which is stored in difG
names(difG)
sapply(difG, names)
ag <- difG
#Extract highly expressed genes from each cell type
ag <- lapply(ag, function(x){
  rownames(x$dg) [x$dg$fc > 1.25 & x$dg$padj < 0.01]
})
sapply(ag, length)
#Take 110 genes for each group (just an arbritary number), which intersect with Class A genes
ag2 <- sapply(ag[c(1,3,2)], function(x){ #Order: 1,3,2 to go in the following order: T cells, B cells, monocytes
  x2 <- intersect(Agen, x)
  x2 <- x2[1:110]
  return(x2)
})
ag2 <- unique(as.vector(ag2))
length(ag2)
head(ag2)

#Order cells by cluster cell type:
plotmap(sN, um=TRUE)
ac <- sN@cpart
#Order of clusters of interest:
ac2 <- c(16,2,17,3,11,18,1,14,20,6,12,10,13,22,4,7,5,9,15,8)
ac3 <- c()
for (i in ac2){
  if (i == ac2[1]) ac3 <- ac[ac == i]
  else ac3 <- c(ac3, ac[ac == i])
}
#ac3 is a vector with the partition from sN@cpart, and with names equal to the cell barcodes

#Plot heatmaps
#Class A genes:
plotmarkergenes(sN,genes=ag2, cells=names(ac3), order.cells=TRUE, noise=FALSE, cluster_rows = FALSE) #expression
plotmarkergenes(sN,genes=ag2, cells=names(ac3), order.cells=TRUE, noise=TRUE, cluster_rows = FALSE) #noise
plotmarkergenes(sN2,genes=ag2,cells=names(ac3), order.cells=TRUE,noise=TRUE,cluster_rows=FALSE) #gene activity

#Class B genes:
#Include hierarchical clustering for these genes
ph <- plotmarkergenes(sN,genes=Bgen, cells=names(ac3), order.cells=TRUE, noise=FALSE) #expression
ph2 <- plotmarkergenes(sN,genes=Bgen[ph$tree_row$order],cells=names(ac3), order.cells=TRUE,noise=TRUE,cluster_rows=FALSE) #noise
ph3 <- plotmarkergenes(sN2,genes=Bgen[ph$tree_row$order],cells=names(ac3), order.cells=TRUE,noise=TRUE,cluster_rows=FALSE) #gene activity

#Class C genes:
#Include hierarchical clustering for these genes (this step can take some minutes)
ph <- plotmarkergenes(sN,genes=Cgen, cells=names(ac3), order.cells=TRUE, noise=FALSE) #expression
ph2 <- plotmarkergenes(sN,genes=Cgen[ph$tree_row$order],cells=names(ac3), order.cells=TRUE,noise=TRUE,cluster_rows=FALSE) #noise
ph3 <- plotmarkergenes(sN2,genes=Cgen[ph$tree_row$order],cells=names(ac3), order.cells=TRUE,noise=TRUE,cluster_rows=FALSE) #gene activity


#Gene set enrichment analysis of marker genes and overlaps with Class A, B or C genes
#Output of differential expression analyses is on difG list
names(difG)
#Extract again marker genes
ag <- difG
ag <- lapply(ag, function(x){
  rownames(x$dg) [x$dg$fc > 1.25 & x$dg$padj < 0.01]
})
sapply(ag, length)
#Subset of genes detected in the matrices of gene activities, gene expression and noise:
k <- intersect(rownames(gene.ac2), rownames(nData))
k <- intersect(k, rownames(eps))
ag2 <- lapply(ag, function(x) intersect(x, k))
sapply(ag2, length)

#Make a list with intersects of marker genes and the groups of class A, B and C
ag3 <- lapply(ag2, function(x){
  x2 <- list()
  x2[[1]] <- x
  x2[[2]] <- intersect(x, Agen) #Overlap: marker genes and Class A genes
  x2[[3]] <- setdiff(x, Agen) #No overlap between marker genes and Class A genes
  x2[[4]] <- intersect(x, Cgen) #Overlap: marker genes and Class C genes
  x2[[5]] <- intersect(x, Bgen) #Overlap: marker genes and Class B genes
  return(x2)
})
sapply(ag3, names)
#Add names to each slot in the ag3 list
x5 <- c("AllMarkers", "OverlapA", "NoOverlapA", "OverlapC", "OverlapB")
ag3 <-lapply(ag3, function(x) {
  names(x) <- x5
  return (x)
})

#Enrichment functions
#Define the universe of genes, as k, the intersect of all detected genes
uni <- k
uni <- bitr(uni, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

#Gene names into ENTREZID
names(ag3)
#T cells:
aT <- lapply(ag3[[1]], function(x){
  bitr(x, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")$ENTREZID
})
sapply(aT, length)
sapply(aT, head)
#Monocytes:
aM <- lapply(ag3[[2]], function(x){
  bitr(x, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")$ENTREZID
})
sapply(aM, length)
sapply(aM, head)
#B cells:
aB <- lapply(ag3[[3]], function(x){
  bitr(x, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")$ENTREZID
})
sapply(aB, length)
sapply(aB, head)

#Enrichment with lists of genes:
xT <- compareCluster(geneClusters=aT, fun = "enrichPathway", organism="human", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, readable=TRUE, universe=uni[,2])
dim(xT)
head(xT)
table(as.data.frame(xT)$Cluster)

xM <- compareCluster(geneClusters=aM, fun = "enrichPathway", organism="human", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, readable=TRUE, universe=uni[,2])
dim(xM)
table(as.data.frame(xM)$Cluster)

xB <- compareCluster(geneClusters=aB, fun = "enrichPathway", organism="human", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, readable=TRUE, universe=uni[,2])
dim(xB)
table(as.data.frame(xB)$Cluster)
#Plots:
dotplot(xT, showCategory=5)
dotplot(xM, showCategory=5)
dotplot(xB, showCategory=5)


##
#Motif enrichment with RcisTarget
#Take the information from ag3 list, computed in the previous section
names(ag3)
sapply(ag3, names)
#Take intersect of marker genes and Class A and B
geneLists <- ag3[[1]][c(2,5)]
geneLists[3:4] <- ag3[[2]][c(2,5)]
geneLists[5:6] <- ag3[[3]][c(2,5)]
x5 <- c("OverA", "OverB")
names(geneLists) <- paste(rep(c("Tc", "Mono", "Bc"), each=length(x5)), x5, sep="_")

# Select motif database to use (i.e. organism and distance around TSS)
data(motifAnnotations_hgnc)
#Motif rankings
#.feather file downloaded from: https://resources.aertslab.org/cistarget/
motifRankings <- importRankings("/path_to_repo/RcisTarget/hg19-tss-centered-10kb-10species.mc9nr.feather")

# Motif enrichment analysis:
motifs <- cisTarget(geneLists, motifRankings, motifAnnot=motifAnnotations_hgnc)
head(motifs)
#Extract relevant information for a couple of columns
motifs2 <- motifs[,c(1:5)]
head(motifs2)
#Extract names of the TF_highConf column
motifs2$TF <- motifs2$TF_highConf %>% gsub(" \\(.+", "", .)

#There are rows with more than one potential TF, separated by ";", separate them
#and filter out genes that are really expressed in the dataset (contained in sN@genes)

#Split rows with double TF's into a separate list
xL <- motifs2$TF %>% strsplit(., "; ")
#Vector with reference genes:
genes_ref <- sN@genes

#Filter slots with at least one TF found in genes_ref
f <- sapply(xL, function(y){
  y2 <- sum(y %in% genes_ref) >0
  return(y2)
})

motifs2 <- motifs2[f,]
df <- motifs2 %>% group_by(geneSet) %>% slice_max(n = 3, order_by = NES) %>% as.data.frame()
df
#There is one duplicated transcription factor, check if both are expressed or not
grep("RELA", genes_ref, value=TRUE)
grep("RELB", genes_ref, value=TRUE)

#Change the duplicated TF
df$TF[2]<- "RELB"

df$geneSet <- factor(df$geneSet, levels=names(geneLists))
x <- unique(df$TF)
df$TF <- factor(df$TF, levels=x)
df$rank <- rep(1:3, times=6)

ggplot(df, aes(rank,TF, fill=NES)) +
  geom_tile() +
  facet_grid(.~geneSet) +
  theme_bw()


#Analysis of signal at the peak levels
DefaultAssay(pbmc) <- "ATAC"
#Make a copy of Seurat objects, in orders to compute different links
pbmc2 <- pbmc
pbmc3 <- pbmc
#Compute links between expression and chromatin accessibility (this step can take a long time)
pbmc2 <- LinkPeaks(
  object = pbmc2,
  peak.assay = "ATAC",
  expression.assay = "SCT"
)

#Create a "Noise" assay into the pbmc3 Seurat object and compute link between noise and chromatin accessibility
#LinkPeaks() function can take a while again
eps <- noise$epsilon
eps[is.na(eps)] <- 1e-05
pbmc3[['Noise']] <- CreateAssayObject(counts = eps)
pbmc3 <- LinkPeaks(
  object = pbmc3,
  peak.assay = "ATAC",
  expression.assay = "Noise",
  gene.use = rownames(eps)
)

#Some exploration:
#Expression - accessibility links
Links(pbmc2)
Links(pbmc2)$zscore %>% summary()
#Noise - accessibility links
Links(pbmc3)
Links(pbmc3)$zscore %>% summary()


#Compute differential accessibility test for specific cell populations
daPL <-list()
daPL$TCells <- FindMarkers( object = pbmc, ident.1 = c("CD4 Naive", "CD4 TCM", "CD8 Naive", "Treg", "CD8 TEM_1",
                                                       "CD8 TEM_2", "CD4 TEM", "gdT"), min.pct = 0.2,
                            test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
daPL$BCells <- FindMarkers( object = pbmc, ident.1 = c("Intermediate B", "Memory B", "Naive B"), min.pct = 0.2,
                            test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
daPL$Mono <- FindMarkers( object = pbmc, ident.1 = c("CD14 Mono", "CD16 Mono"), min.pct = 0.2,
                          test.use = 'LR', latent.vars = 'atac_peak_region_fragments')

#Extract peaks with increased / decreased accessibility
#Select a slot of the daPL list
names(daPL)
xL <-daPL$TCells
#Increased accessibility
x <- rownames (xL) [xL$avg_log2FC > log2(1.25) & xL$p_val_adj < 0.001]
daTC <- StringToGRanges(x, sep = c("-", "-"))
daTC$Accessibility <- "Open"
#Decreased accessibility
x <- rownames (xL) [xL$avg_log2FC < log2(0.8) & xL$p_val_adj < 0.001]
daTC2 <- StringToGRanges(x, sep = c("-", "-"))
daTC2$Accessibility <- "Closed"
#Combine both sets of peaks
daTC <- c(daTC, daTC2)
rm(daTC2)

#Increased/decreased accessibility in Monocytes
xL <-daPL$Mono
x <- rownames (xL) [xL$avg_log2FC > log2(1.25) & xL$p_val_adj < 0.001]
daMC <- StringToGRanges(x, sep = c("-", "-"))
daMC$Accessibility <- "Open"

x <- rownames (xL) [xL$avg_log2FC < log2(0.8) & xL$p_val_adj < 0.001]
daMC2 <- StringToGRanges(x, sep = c("-", "-"))
daTM2$Accessibility <- "Closed"
daMC <- c(daMC, daMC2)
rm(daMC2)

#Make plots of genomic tracks
#Chose some genes of interest
g <- "CD28"
Links(pbmc2)[Links(pbmc2)$gene == g]
Links(pbmc3)[Links(pbmc3)$gene == g]

#Check of the whole coverage of links:
#Links: expression - chrom. accessibility
CoveragePlot(pbmc2, region = g, features = g, assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE,
             idents = c("CD4 Naive", "CD8 Naive", "CD14 Mono", "Memory B"),
             extend.upstream = 5e+05, extend.downstream = 5e+05,
             ranges = daTC,
             ranges.title = "DA test",
             ranges.group.by="Accessibility"
)
#Links: noise - chrom. accessibility
CoveragePlot(pbmc3, region = g, features = g, assay = 'ATAC', expression.assay = 'Noise', peaks = TRUE,
             idents = c("CD4 Naive", "CD8 Naive", "CD14 Mono", "Memory B"),
             extend.upstream = 5e+05, extend.downstream = 5e+05,
             ranges = daTC,
             ranges.title = "DA test",
             ranges.group.by="Accessibility"
)

#Adjust the genomic region to be closer to the gene of interest:
#Links: expression - chrom. accessibility
CoveragePlot(pbmc2, region = g, features = g, assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE,
             idents = c("CD4 Naive", "CD8 Naive", "CD14 Mono", "Memory B"),
             extend.upstream = 15000, extend.downstream = 25000,
             ranges = daTC,
             ranges.title = "Open peaks",
             ranges.group.by="Accessibility"
)
#Links: noise - chrom. accessibility
CoveragePlot(pbmc3, region = g, features = g, assay = 'ATAC', expression.assay = 'Noise', peaks = TRUE,
             idents = c("CD4 Naive", "CD8 Naive", "CD14 Mono", "Memory B"),
             extend.upstream = 15000, extend.downstream = 25000,
             ranges = daTC,
             ranges.title = "DA test",
             ranges.group.by="Accessibility"
)


#Chose another gene
g <- "AKAP13"
Links(pbmc2)[Links(pbmc2)$gene == g]
Links(pbmc3)[Links(pbmc3)$gene == g]

#Links: expression - chrom. accessibility
CoveragePlot(pbmc2, region = g, features = g, assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE,
             idents = c("CD4 Naive", "CD8 Naive", "CD14 Mono", "Memory B"),
             extend.upstream = 15000, extend.downstream = 25000,
             ranges = daMC,
             ranges.title = "DA test",
             ranges.group.by="Accessibility"
)
#Links: noise - chrom. accessibility
CoveragePlot(pbmc3, region = g, features = g, assay = 'ATAC', expression.assay = 'Noise', peaks = TRUE,
             idents = c("CD4 Naive", "CD8 Naive", "CD14 Mono", "Memory B"),
             extend.upstream = 15000, extend.downstream = 25000,
             ranges = daMC,
             ranges.title = "DA test",
             ranges.group.by="Accessibility"
)

