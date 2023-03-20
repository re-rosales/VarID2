## VarID2 is part of the RaceID (>= v0.2.4) package.
## This script contains the code to reproduce the initial analysis and benchmarking of VarID2 method,
## contained in Figures 1 and S1 of Rosales-Alvarez et al., manuscript
## Data from Tusi, et al. Population snapshots predict early haematopoietic and erythroid hierarchies. Nature 555, 54–60 (2018).

## The analysis in this script was done with RaceID v0.2.6. For further versions, cluster numbers and parameter settings potentially change.

# Setup python environment in order to use leiden clustering
# It requires to install reticulate R package and the python "leidenalg" and "igraph" modules in the indicated python environment
reticulate::use_python("./anaconda3/bin/python", required = TRUE)
# Corroborate "leidenalg" and "igraph" are available, as well as python
reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()

# Load required libraries
library(RaceID)
library(Matrix)
library(MASS)
library(FNN)
library(pheatmap)
library(ggplot2)
library(scales)
source("VarID2_Fig1_SideFunctions.R")

# Load input data (Tusi et al., 2018)
#This file can be downloaded from: https://github.com/dgrun/VarID_analysis/blob/master/inputData_hematopoiesis.rds
inputData <- readRDS("inputData_hematopoiesis.rds")

# General pipeline to create a sC object, quality filtering and clustering by building a knn pruned network.
sN <- SCseq(inputData$wt)
sN <- filterdata(sN,mintotal=1000,CGenes=grep("^(mt|Rp(l|s)|Gm\\d)",rownames(inputData$wt),value=TRUE))

expData  <- as.matrix(getExpData(sN))
res      <- pruneKnn(expData,no_cores=5,knn=100)
cl    <- graphCluster(res,pvalue=0.01,use.leiden=TRUE,leiden.resolution=1)
probs <- transitionProbs(res,cl)

sN <- updateSC(sN,res=res,cl=cl)
sN <- comptsne(sN)
sN <- compumap(sN)
plotmap(sN, um=TRUE)


## Exploring components of variability: coefficient of variation as a function of the mean
#Retrieve the matrix of prunned knn with indexes
pnn <- t(res$NN)[, -1] * (t(res$pvM) > 0.01)
pnn <- cbind(t(res$NN)[, 1], pnn)
#Select a neighborhood with a central cell i
#Neighborhood 1, selected in the manuscript:
i <- "bBM17"
f <- pnn[i,]
f <- f[f !=0]

#Obtain mean UMI counts per gene and variance
expData  <- as.matrix(getExpData(sN)) #extract matrix with raw UMI counts
m  <- apply(expData[,f], 1, mean)
v  <- apply(expData[,f], 1, var)

#Fit a gamma distribution for the total UMI counts per cell within the neighborhood with central cell i
k <- colSums(expData[,f])
gfit <- tryCatch(fitdistr(k/mean(k), "gamma"), error = function(err) return(NA))
#Dispersion parameter rt:
rt <-  gfit$estimate["shape"]

#Get mean and coefficient of variation into logarithmic scale
l_m <- log2(m)
l_cv <- 0.5*log2(v) - log2(m)
plot(l_m, l_cv, col="gray", pch=20, xlab="log2 mean expression", ylab="log2 CV")
#Poissonnian noise
abline(0, -0.5, lwd=2) 
#Total UMI count variability, corresponding to the variance of the gamma fit
gV <- gfit$estimate["shape"]/gfit$estimate["rate"]^2
abline(h = 0.5 * log2(gV), col = "orange", lwd=2)

#Compute explained variance and add it into the dot plot
xc <- m[order(m, decreasing = FALSE)]
l_xc <- log2(xc)
cv.explained <- sigmaNB(xc, rt)/xc
l_cv.explained <- log2(cv.explained)

lines(l_xc, l_cv.explained, col = "purple", lwd=2)
legend("topright", legend=c("Poissonian noise", "Total UMI count var.", "Explained var."), col=c("black", "orange","purple"),
       lty=1, lwd=2, bty="n", main=paste("Neighborhood", i, sep=" #"))

#Repeat the same code for another neighborhood:
#Neighborhood 2, selected in the manuscript:
j <- "bBM46"
f <- pnn[j,]
f <- f[f !=0]
#Mean and variance:
m2  <- apply(expData[,f], 1, mean)
v2  <- apply(expData[,f], 1, var)
#Gamma fit of UMIs counts variance and rt dispersion parameter inference
k2 <- colSums(expData[,f])
gfit2 <- tryCatch(fitdistr(k2/mean(k2), "gamma"), error = function(err) return(NA))
rt2 <-  gfit2$estimate["shape"]

#Get mean and coefficient of variation into logarithmic scale
l_m2 <- log2(m2)
l_cv2 <- 0.5*log2(v2) - log2(m2)
plot(l_m2, l_cv2, col="gray", pch=20, xlab="log2 mean expression", ylab="log2 CV")
#Poissonnian noise
abline(0, -0.5, lwd=2) 
#Total UMI count variability, corresponding to the variance of the gamma fit
gV2 <- gfit2$estimate["shape"]/gfit2$estimate["rate"]^2
abline(h = 0.5 * log2(gV2), col = "orange", lwd=2)
#Explained variance:
xc2 <- m2[order(m2, decreasing = FALSE)]
l_xc2 <- log2(xc2)
cv.explained2 <- sigmaNB(x2c, rt2)/xc2
l_cv.explained2 <- log2(cv.explained2)

lines(l_xc2, l_cv.explained2, col = "purple", lwd=2)
legend("topright", legend=c("Poissonian noise", "Total UMI count var.", "Explained var."), col=c("black", "orange","purple"),
       lty=1, lwd=2, bty="n", main=paste("Neighborhood", j, sep=" #"))


#Comparing number of UMI counts and detected features per barcode for each cellular neighborhood
#Neighborhood 1:
i <- "bBM17"
f <- pnn[i,]
f <- f[f !=0]
#UMI counts per barcode:
um1 <- colSums(expData[,f])
#Number of features with at least 1 UMI count in each barcode:
nf1 <- apply(expData[,f], 2, function(x) sum(x>0))
#Neighborhood 2:
j <- "bBM46"
f <- pnn[j,]
f <- f[f !=0]
um2 <- colSums(expData[,f])
nf2 <- apply(expData[,f], 2, function(x) sum(x>0))

#Make a data fame to make violin plots
df <- data.frame(UMIs=um1, nFeat=nf1, neigborhood="Neig.1")
df <- rbind(df, data.frame(UMIs=um2, nFeat=nf2, neigborhood="Neig.2"))
#Order the names of the neighborhoods for better visualization
df$neigborhood <- factor(df$neigborhood, levels=c("Neig.2", "Neig.1"))

#Plot number of UMI counts per barcodes
ggplot(df, aes(x=neigborhood, y=UMIs, fill=neigborhood)) +
  scale_fill_manual(values=c("#ff7f7f", "#7f7fff")) +
  geom_violin(scale="width", trim=FALSE) +
  coord_flip() +
  theme_pubr() +
  ylim(0,max(df$UMIs)) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black", size=2) +
  ggtitle("UMI counts") +
  theme(legend.position = "none")

#Plot number of detected features per barcodes
ggplot(df, aes(x=neigborhood, y=nFeat, fill=neigborhood)) +
  scale_fill_manual(values=c("#ff7f7f", "#7f7fff")) +
  geom_violin(scale="width", trim=FALSE) +
  coord_flip() +
  theme_pubr() +
  ylim(0,max(df$nFeat)) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black", size=2) +
  ggtitle("Detected features") +
  ylab("Number of detected features") +
  theme(legend.position = "none")



## Gamma fit of local neighborhoods
#Compute the normalization parameter beta for each cellular neighborhood.
#Beta parameter is defined as the total UMI counts per cell within a neighborhood divided by the mean UMI counts
#Neighborhood 1:
i <- "bBM17"
f <- pnn[i,]
f <- f[f !=0]
k1 <- colSums(expData[,f])
x1 <- k1/mean(k1)

#Fit a gamma distribution
gfit1 <- tryCatch(fitdistr(x1, "gamma"), error = function(err) return(NA))
#Dispersion parameter rt:
rt1 <-  gfit1$estimate["shape"]
#Density function:
xr1 <- seq(0,max(x1),length.out=100)
d1 <- dgamma(xr1,shape=rt1,rate=rt1)

#Neighborhood 2:
j <- "bBM46"
f <- pnn[j,]
f <- f[f !=0]
k2 <- colSums(expData[,f])
x2 <- k2/mean(k2)

#Fit a gamma distribution
gfit2 <- tryCatch(fitdistr(x2, "gamma"), error = function(err) return(NA))
#Dispersion parameter rt:
rt2 <-  gfit1$estimate["shape"]
#Density function:
xr2 <- seq(0,max(x2),length.out=100)
d2 <- dgamma(xr2,shape=rt2,rate=rt2)

#Histogram of calculated beta values and gamma fit
h <- hist(c(x1, x2),breaks=30,plot=FALSE)
h1 <- hist(x1,breaks=h$breaks,plot=FALSE)
h2 <- hist(x2,breaks=h$breaks,plot=FALSE)

plot(h1,border="white",col=alpha("blue", .4),freq=FALSE, xlim=range(c(x1,x2)), ylim=c(0,max(c(h1$density, h2$density))) )
plot(h2,border="white",col=alpha("red", .4),freq=FALSE, add=TRUE)
lines(xr1, d1, col="blue")
lines(xr2, d2, col="red")


#Estimate rt (technical dispersion parameter), as long as individual gene noise across all cell neighborhoods in the hematopoietic dataset
#Use compTBNoise function, with an arbitrary gamma value for now.
nres <- compTBNoise(res,expData,pvalue=.01,no_cores=5,gamma=2,x0=0,lower=0,upper=100)
#Retrieve rt parameter from the output list
rtH <- nres$rt
plotfeatmap(sN, rtH, um=TRUE, cex=0.5)

## Cauchy distribution as a prior distribution for the residual biological noise
#Location parameter is set to 0
#Exploration of different gamma values (scale parameter):
x <- seq(from=0.1, to=10, length.out=100)
plot(x, dcauchy(x, location = 0, scale = 5), type="l", ylim=c(0,1.5), lwd=2, col="cyan")
lines(x, dcauchy(x, location = 0, scale = 2), lwd=2, col="pink")
lines(x, dcauchy(x, location = 0, scale = 1), lwd=2, col="blue")
lines(x, dcauchy(x, location = 0, scale = 0.5), lwd=2, col="#A10505")
lines(x, dcauchy(x, location = 0, scale = 0.2), lwd=2, col="#ff7f24")
abline(h=0, lty=5, lwd=2, col="gray")

cols <- c("cyan","pink", "blue","#A10505", "#ff7f24", "gray")
y <- c(paste("Gamma", c(5,2,1, 0.5,0.2), sep=":"), "Horizontal line")
y2 <- c(rep(1,5), 5)
legend("topright", legend=y, lty=y2, col=cols, bty="n")



## Data simulation
#Take as reference estimates from the hematopoietic dataset.
#Get raw UMI counts
expData <- getExpData(sN)
#Compute mean, variance and dispersion parameter per gene
m  <- apply(expData, 1, mean)
v  <- apply(expData, 1, var)
r  <- 1/(m**2/(v - m))
r  <- r[!is.na(r) & r > 0]

#Matrix of 100 nearest neighbors:
k  <- t(res$NN)
#Sample 200 cells
n  <- 200
set.seed(123)
sa <- sample(1:ncol(expData),n)
#Obtain mean and variance across the 200 random neighborhoods
for ( i in 1:length(sa) ){
  n <- k[sa[i],]
  ml <- if ( i == 1 ) apply(expData[,n],1,mean) else cbind(ml,apply(expData[,n],1,mean))
  vl <- if ( i == 1 ) apply(expData[,n],1,var) else cbind(vl,apply(expData[,n],1,var))
}
ml <- cbind(ml,apply(expData,1,mean))
vl <- cbind(vl,apply(expData,1,var))
#Estimate dispersion parameter
rl  <- 1/(ml**2/(vl - ml))
rl  <- rl[!is.na(rl) & rl > 0]

#Generate a vector with mean expression, by sampling values from m vector as a reference
set.seed(123)
nn <- sample(m,length(m)*10,replace=TRUE)

#Reference vector of values for the dispersion parameter r
#Define three levels of dispersion, based on quantiles from rl vector.
#Quantiles: 0.8, 0.4 and 0.2, for low, medium and high dispersion levels
rs <- c(
  rep(quantile(rl,.8),round(length(nn)*.8,0)),
  rep(quantile(rl,.4),round(length(nn)*.1,0)),
  rep(quantile(rl,.1),length(nn) - round(length(nn)*.8,0) - round(length(nn)*.1,0))
)
rs <- round(rs,2)

#Generate alpha values, corresponding to per-cell variability in UMI counts
par <- 2 
nb <- 100
alpha <- rgamma(nb,shape=par,rate=par)

#Data simulation, taking as reference the previous estimates: mean, dispersion parameter r and per-cell variability in UMI counts; contained in nn, rs and alpha, respectively
d  <- apply(cbind(nn,rs),1,function(x){ sapply(alpha,function(y) rnbinom(1,size=x[2],mu=y*x[1]))})
d <- t(d)
rownames(d) <- paste("Gene",rs,1:nrow(d),sep="_")
colnames(d) <- paste("Cell",1:ncol(d),sep="_")
Z  <- d

#Visualization of the simulated data
#Plot variance as a function of the mean, highlighting different levels of noise
M <- apply(Z,1,mean)
V <- apply(Z,1,var)
lev <- rev(sort(unique(rs)))
f1 <- rs == lev[1]
f2 <- rs == lev[2]
f3 <- rs == lev[3]
f <- M > 0
plot(log2(M)[f],log2(V)[f],type="n", xlab="Mean expression", ylab="Variance")
points(log2(M)[f & f1],log2(V)[f & f1],cex=.5,col="red")
points(log2(M)[f & f2],log2(V)[f & f2],cex=.5,col="turquoise")
points(log2(M)[f & f3],log2(V)[f & f3],cex=.5,col="blue")


#Tests for different gamma values (for the Cauchy prior)
#List for storing the main results, after optimization:
W    <- list()
#Gamma values to test:
#1000 gives a flat prior, so it is like performing maximum likelihood estimation
GAMMA <- c(0.2,0.5,1,2,5,1000)
#Use fitNBtb function for each gamma value. This step can take a some minutes, depending on the computational resources
for (gamma in GAMMA){
  x <- fitNBtb(Z,gamma,x0,lower=0, upper=100)
  W[[paste("gamma",gamma,sep=":")]] <- x
}

#Extract the gene expression noise values, contained in the "epsilon" slot
VFC <- lapply(W, function(x) x$epsilon)
#VFC  <- list() #It contains only epsilon estimates
#X0   <- 0:5/2
#names(W)
#names(W$"x0:0")
#sapply(W$"x0:0", dim)
#sapply(W$"x0:0", head)
#length(rs)
#names(VFC)

#Restrict quality measurements to gene with mean expression higher than one
M <- apply(Z,1,mean)
f  <- M > 1
#Create matrices to store some metrics:
lev <- rev(sort(unique(rs)))
MED <- matrix(NA, nrow=length(lev), ncol=length(GAMMA))
rownames(MED) <- paste("eps", 1/lev, sep=":")
colnames(MED) <- colnames(GAMMA)
DEV <- DELTA <- MED

for (gamma in GAMMA){
  vfc <- VFC2[[paste("gamma",gamma,sep=":")]]
  vfc <- vfc + .01
  for ( i in 1:length(lev) ){
    g <-  rs == lev[i]
    MED[i,gamma]   <- median(vfc[g & f],na.rm=TRUE)/(1/lev[i])  #Ratio: median of epsilon / Ground truth epsilon (inverse of lev values)
    DELTA[i,gamma] <- sqrt(var(log(vfc)[!is.na(f) & f & g]))  #standard deviation of epsilon
    DEV[i,gamma] <- sum( log(vfc)[g] > log(1/lev[i]) + delta[i], na.rm=TRUE) #Number of outliers
    }
}

#Plot: ratio between median of epsilon / Ground truth epsilon
pheatmap(MED2,cluster_rows=FALSE,cluster_cols=FALSE, fontsize=15)
#Plot: standard deviation of epsilon
pheatmap(DELTA2,cluster_rows=FALSE,cluster_cols=FALSE, fontsize=15)
#Plot: number of outliers
pheatmap(DEV2,cluster_rows=FALSE,cluster_cols=FALSE, fontsize=15)


#Comparison of maximum a posteriori (MAP) and maximum likelihood (ML) estimates
#MAP estimates (with gamma = 1)
w   <- W[["gamma:1"]]
vfc <- VFC[["gamma:1"]]
#ML estimates
wm   <- W[["gamma:1000"]]
vmc <- VFC[["gamma:1000"]]

lev <- rev(sort(unique(rs)))
f1 <- rs == lev[1]
f2 <- rs == lev[2]
f3 <- rs == lev[3]
f <- w$mu > 0

#Setup common y limit for better visualization
yl <- range(c(log2(vfc)[f & log2(vfc) > -Inf], log2(vmc)[f & log2(vmc) > -Inf]))
#Plot genes with high noise levels, whow index is indicated in f3 vector:
plot(log2(wm$mu)[f & f3],log2(vmc)[f & f3],pch=20,col="black",xlab="Mean expression",ylab="Epsilon", ylim=a)
points(log2(w$mu)[f & f3],log2(vfc)[f & f3],cex=1.5,col="red")
legend("bottomleft", c("MLE", "MAP"), pch=c(20,1), col=c("black", "red"), bty="n")

#Plot epsilon estimates as a function of the mean
#Segregate estimates based on the noise level: low, medium and high, indicated in f1, f2 and f3 vectors
plot(log2(w$mu)[f],log2(vfc)[f],type="n",xlab="log2 mean expression",ylab="log2 epsilon")
points(log2(w$mu)[f & f1],log2(vfc)[f & f1],cex=.5,col="red")
points(log2(w$mu)[f & f2],log2(vfc)[f & f2],cex=.5,col="turquoise")
points(log2(w$mu)[f & f3],log2(vfc)[f & f3],cex=.5,col="blue")
#Ground truth epsilon estimates, as the inverse of values contained in lev vector:
gt_eps <- round(1/lev, 2)
abline(h=log2(gt_eps[1]),lt=2,col="red")
abline(h=log2(gt_eps[2]),lt=2,col="turquoise")
abline(h=log2(gt_eps[3]),lt=2,col="blue")
#Median values of the epsilon estimates:
abline(h=log2(median(vfc[f & f1],na.rm=TRUE)),col="red")
abline(h=log2(median(vfc[f & f2],na.rm=TRUE)),col="turquoise")
abline(h=log2(median(vfc[f & f3],na.rm=TRUE)),col="blue")

legend("bottomleft", legend=c("Ground truth", paste("e =", gt_eps[3], sep=" "), paste("e =", gt_eps[2], sep=" "),
                              paste("e =", gt_eps[1], sep=" ")), col=c("white", "blue", "turquoise", "red"),
       lty=2, bty="n")
legend("bottomright", legend=c("epsilon estimates, \n median values", "high", "medium", "low"),
       col=c("white","blue", "turquoise", "red"),
       lty=1, bty="n")


#Explore a quality control for genes
#Set a threshold of values for mean expression:
thr <- c(0,0.1,0.2,0.3,0.4,0.5,1)
#Use values from w data frame, containing MAP estimates with gamma=1
head(w)
ng <- list(low=w[f1,], medium=w[f2,], high=w[f3,])
#Add ground truth noise levels contained in gt_eps vector
for (i in 1:length(ng)){
  ng[[i]]$noiseGT <- gt_eps[i]
}

#Make a test, and quantify two measurements:
#Fraction of genes passing the mean threshold values
#Percentage of genes within 1-fold change interval around the ground truth noise level
ngf <- lapply(ng, function(x){
  muX <- x$mu
  epsilonX <- x$epsilon
  fc <- unique(x$noiseGT)
  #Upper and lower threshold values from the ground truth:
  fc1 <- fc*2
  fc2 <- fc*0.5
  df <- matrix(NA, nrow=2, ncol=length(thr))
  for (i in 1:length(thr)){
    f <- muX >= thr[i]
    #Percentage of genes over the mean threshold:
    df[1,i] <- round((sum(f)/length(muX)) *100, 2) 
    #Percentage of genes within the 1-fold change interval
    df[2,i] <- round((sum(epsilonX[f] < fc1 & epsilonX[f] > fc2) / sum(f))*100, 2)
  }
  colnames(df) <- paste("mean >=", thr)
  rownames(df) <- c("GenesOverExprThres", "GenesWithinInterval")
  return(df)
})

#Arrange the quantifications based on the measurement, in order to visualize the scores across the three groups of genes
ng2f <- list(GenesOverExprThres = t(sapply(ngf, function(x) x[1,])),
             GenesWithinInterval = t(sapply(ngf, function(x) x[2,])) )
pheatmap(ng2f[[1]],cluster_rows=FALSE,cluster_cols=FALSE, main="Fraction Genes > mean expr.")
pheatmap(ng2f[[2]],cluster_rows=FALSE,cluster_cols=FALSE, main="Fraction Genes < 2FC")


#Comparison with BASiCS method
library(BASiCS)
#Use the simulated data contained in Z matrix as an input
Counts <- as.matrix(Z)
Data <- SingleCellExperiment(list(counts=Counts))
#Assign cells to two different batches:
x <- ceiling(ncol(Counts) / 2)
colData(Data)$BatchInfo<-c(rep(1:2, times=a)[1:ncol(Counts)])
#Next command takes a long time, depending on the computational resources
Chain <- BASiCS_MCMC(Data = Data,  N = 20000, Thin = 20, Burn = 10000,PrintProgress = FALSE, Regression = TRUE, WithSpikes = FALSE)
#Extract the summaries of the posterior distributions and the median values from epsilon parameter
ChainSummary <- Summary(Chain)
eps <- displaySummaryBASiCS(ChainSummary, Param = "epsilon")
eps <- eps[,1]

#Add a pseudocount to MAP and ML estimates, contained in vfc and vmc, respectively
y <- vfc + .01
ym <- vmc + .01
#Common y limits for better visualization on the plots:
yl <- range(c(log2(y), log2(ym)) )

#Comparison BASiCS and MAP estimates
plot(eps,log2(y),type="n", main="Comparison MAP and BASICS", xlab="e BASiCS", ylab="log2 (e MAP + 0.01)", ylim=yl)
points(eps[f1],log2(y[f1]),col="red",cex=.5)
points(eps[f2],log2(y[f2]),col="turquoise",cex=.5)
points(eps[f3],log2(y[f3]),col="blue",cex=.5)
legend("bottomright", legend=c("Noise levels", "High", "Medium", "Low"),
       col=c("white", "blue", "turquoise", "red"),
       pch=1, bty="n")
cor(y, eps, method="pearson", use="complete.obs")

#Comparison BASiCS and MLE estimates
plot(eps,log2(ym),type="n", main="Comparison MLE and BASICS", xlab="e BASiCS", ylab="log2 (e MLE + 0.01)", ylim=yl)
points(eps[f1],log2(ym[f1]),col="red",cex=.5)
points(eps[f2],log2(ym[f2]),col="turquoise",cex=.5)
points(eps[f3],log2(ym[f3]),col="blue",cex=.5)
legend("bottomright", legend=c("Noise levels", "High", "Medium", "Low"),
       col=c("white", "blue", "turquoise", "red"),
       pch=1, bty="n")
cor(ym, eps, method="pearson", use="complete.obs")


#Comparison of individual (1D) and simultaneous (2D) parameter estimation
#Individual estimation (1D): inference of epsilon (noise) estimates, while mean expression is directly calculated from the dataset. It is the conventional VarID2 approach
#Simultaneous estimation (2D): inference of both epsilon and mean expression

#Use data from hematopoetic cells (Tusi et al., 2018)
#Use the matrix of k-nearest neighbors
k <- t(res$NN)
#Select a radom neighborhood
j <- 100
expData <- getExpData(sN)
Z <- expData[,k[j,]]

gamma  <- 1 
x0    <- 0
upper <- 100
lower <- 0

#Individual estimation (1D)
w1 <- fitNBtbTest(Z,gamma,x0,lower,upper)
#Simultaneous estimation (2D)
w2 <- fitNBtbTest(Z,gamma,x0,lower,upper,grad=FALSE)

#Comparison of noise values:
plot(w1$epsilon + 1e-4,w2$epsilon + 1e-4,log="xy", xlab="epsilon with calculated mu", ylab="epsilon with inferred mu",pch=20,col="grey")
abline(0,1)
f <- !is.na(w2$epsilon) & (w2$epsilon > 2* w1$epsilon | w2$epsilon < 0.5* w1$epsilon)
points(w1$epsilon[f] + 1e-4,w2$epsilon[f] + 1e-4,col="red")
legend("bottomright", paste(sum(f), " out of ", nrow(w1)," outliers (",round(sum(f)/nrow(w1)*100,1),"%)",sep=""))

#Comparison of mean expression values:
plot(w1$mu + 1e-2,w2$mu + 1e-2, log="xy", xlab="calculated mu", ylab="inferred mu",pch=20,col="grey")
points(w1$mu[f] + 1e-2,w2$mu[f] + 1e-2,col="red")
abline(0,1)
