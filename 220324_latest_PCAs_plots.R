#########################
##This script was written by Luke Archer (2024) to produce PCA plots in R.
## Load vcf, run PCA, calculate K-means clustering, calculate matrix of geneic distances, run AMOVA
## Filip Kolar 2017, further edits by Sian Bray 2018 and Levi Yant 2022,3

setwd("/Users/lukearcher/Desktop/")

options(warn=1)

library(adegenet)
library(adegraphics) #not strictly necessary for all of this (hombrew r installs will interfere)
library(vcfR)
library(pegas)
library(StAMPP)
library(ade4)
library(MASS)

##################### 
# MODIFIED FUNCTIONS for TETRAPLOIDS

# This block is a function for conversion from vcfR object to genlight in tetraploids and hexaploids: note this changed function is not necessary for LIFE4141 assignment, but it may be helpful later.
vcfR2genlight.tetra <- function (x, n.cores = 1) 
{
  bi <- is.biallelic(x)
  if (sum(!bi) > 0) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
    warning(msg)
    x <- x[bi, ]
  }
  x <- addID(x)
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0
  x[x == "0|1"] <- 1
  x[x == "1|0"] <- 1
  x[x == "1|1"] <- 2
  x[x == "0/0"] <- 0
  x[x == "0/1"] <- 1
  x[x == "1/0"] <- 1
  x[x == "1/1"] <- 2
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/1"] <- 1
  x[x == "0/0/0/0/1/1"] <- 2
  x[x == "0/0/0/1/1/1"] <- 3
  x[x == "0/0/1/1/1/1"] <- 4
  x[x == "0/1/1/1/1/1"] <- 5
  x[x == "1/1/1/1/1/1"] <- 6
  if (requireNamespace("adegenet")) {
    x <- new("genlight", t(x), n.cores = n.cores)
  }
  else {
    warning("adegenet not installed")
  }
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  return(x)
}



## -------------------------------------------------------
### This is a patch for MUCH MUCH faster PCA calculation on genlight objects
# see https://github.com/thibautjombart/adegenet/pull/150
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks..
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar
  # to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS
  # to support thousands of samples, this could be
  # replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  ## scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}

# ---------------------------------------------------------

#############################
# IMPORT SNP data from VCF
#import the SNP data for lyrata VCF
lyrata_vcf <- read.vcfR("220324_tetraploids_only.vcf.gz")

###convert the lyrata VCF to genlight object
#convert to genlight 	
aa.genlight <- vcfR2genlight.tetra(lyrata_vcf)                           ## use the modified function vcfR2genlight.tetra at the end of the file
locNames(aa.genlight) <- paste(lyrata_vcf@fix[,1],lyrata_vcf@fix[,2],sep="_")   # add real SNP.names
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3)

#check ========IMPORTANT TO CHECK========
aa.genlight
indNames(aa.genlight)
ploidy(aa.genlight)
####################################

################################
#   PCA 
##################################
# run a PCA
pca.1 <- glPcaFast(aa.genlight, nf=40) # use the modified function glPcaFast at the end of the file

glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks..
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar
  # to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS
  # to support thousands of samples, this could be
  # replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  ## scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}
######
###PLOTTING THE LATEST DATA
scatter(pca.1, posi='bottomleft')
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis

col <- funky(10)
s.class(pca.1$scores, pop(aa.genlight), xax=1, yax=2, col=transp(col,.6), 
        ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F)

##plot the first 2 principal components; colour by population
plot(pca.1$scores[,1], pca.1$scores[,2], 
     cex=2, pch=20, bg="beige",col=aa.genlight$pop, 
     xlab="PC1: 13.55% variance", 
     ylab="PC2: 5.91% variance", 
     main="PCA on Arabidopsis lyrata tetraploids (57,613 SNPs)")
legend("topright", 
       legend=unique(aa.genlight$pop), 
       pch=20, 
       col=funky(10))

########Phylogenetic tree building with Ape######################
##make a phylogenetic tree of the data to see how the samples group together
tree1<-nj(dist(as.matrix(aa.genlight)))
tree1
#plot the tree
plot(tree1,type='cladogram',cex=0.7)
title("Unrooted NJ tree of tetraploid Arabidopsis lyrata")
#######PHYLOGENETIC analysis
plot(tree1, typ="phylogram", show.tip=TRUE,cex=0.8)
tiplabels(pch=20, col=aa.genlight$pop, cex=4)
title("NJ tree of tetraploid A. lyrata")

################################
##make a different PCA plot with colours
myCol <- colorplot(pca.1$scores,pca.1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca.1$eig[1:30],2,1,2, posi="topright", inset=.05, ratio=.3)

######SplitsTREE phylogenetics###############
################################################################################
#===============================================================================
#  distance-based analyses     -------------------------------------------------

# Calculate Nei's distances between individuals/pops
#Ana# Note here raw data is used, not the corrected for NAs
# ---

aa.D.ind <- stamppNeisD(aa.genlight, pop = FALSE) # Nei's 1972 distance between indivs
# export matrix - for SplitsTree
stamppPhylip(aa.D.ind, file="aa.indiv_Neis_distance_4ds.phy.dst")

aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE)   # Nei's 1972 distance between pops
# export matrix - for SplitsTree
stamppPhylip(aa.D.pop, file="aa.pops_Neis_distance_4ds.phy.dst") 

#Ana# the stamppNeisD  is having problem identifying populations here.
# if you use in line 240 pop(aa.genlight) <-substr(indNames(aa.genlight),1,6)
# matrices aa.D.ind and aa.D.pop are the same! what you are calculating is INDIVIDUALS
# since "6" corresponds to full name including the sample nr (1-5)



### create the dist objects                 #Ana# this alters the matrices above
colnames(aa.D.ind) <- rownames(aa.D.ind)
aa.D.ind.dist <-as.dist(aa.D.ind, diag=T)
attr(aa.D.ind.dist, "Labels") <-rownames(aa.D.ind)   # name the rows of a matrix  

colnames(aa.D.pop) <- rownames(aa.D.pop) 
aa.D.pop.dist <-as.dist(aa.D.pop, diag=T)
attr(aa.D.pop.dist, "Labels") <-rownames(aa.D.pop)   # name the rows of a matrix  



#Ana# Table with distances between individuals and populations is now ready!
aa.D.ind # individuals
aa.D.pop # populations
Dgen_ind <- aa.D.ind.dist
Dgen_pop <- aa.D.pop.dist

###################################
##Poly_SFS data
new_lyrata_sfs<-read.csv("220324_SFS.csv",sep=",")
head(new_lyrata_sfs)

new_lyrata_sfs<-c(10479,2333,2432,2348,2336,2289,2210,2186,2088,2004,1798,1746,1673,1439,1361,1208,1037,954,849,709,615,586,509,490,429,391,363,344,359,255,271,251,237,223,227,206,178,191,190,186,171,145,147,139,142,137,125,138,102,111,99,111,113,97,88,92,81,92,75,86,91,67,67,73,63,82,60,61,44,63,62,53,47,52,60,56,43,53,48,55,37,51,36,34,41,35,37,30,38,35,34,35,17,26,27,31,19,30,29,35,19,27,32,21,17,20,28,27,15,20,21,23,19,17,19,15,19,18,14,18,15,23,12,16,8,15,13,13,8,8,14,16,14,15,12,12,11,8,16,13,12,14,13,11,9,12,16,9,15,11,6,13,13,9,4,11,10,6,10,10,4,10,11,12,13,3,9,11,7,14,7,9,14,8,7,12,12,8,8,7,8,7,9,10,4,10,9,8,10,12,7,17,12,16,11,7,9,11,25,35,83)
number_sites<-c(1:201)
prop<-colSums(lyrata_df) ##55519
lyrata_df$proportion<-new_lyrata_sfs/prop
lyrata_df$bins<-1:201

head(lyrata_df)

lyrata_df<-as.data.frame(new_lyrata_sfs,number_sites)
head(lyrata_df)
class(lyrata_df$new_lyrata_sfs)

ggplot(lyrata_df,aes(x=number_sites,y=new_lyrata_sfs))+
  geom_point()+
  theme_bw()+
  labs(title="Arabidopsis lyrata tetraploids site frequency spectrum",
       x="Derived allele frequency",y="Count")+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor.x=element_blank())

################
##pairwise Fst betweeen BZD and OCH 
library(tidyverse)
BZD_OCH_Fst<-read_tsv("BZD_OCH_Fst_output.fst")

head(BZD_OCH_Fst)

ggplot(BZD_OCH_Fst,aes(x=`1317`,y=`0.193548`,colour=threshold))+
  geom_point(alpha=0.4)+
  scale_colour_brewer(palette="Set2")+
  geom_hline(yintercept=0.55,linetype=2,colour=1)+
  theme_bw()+
  labs(title="Fst scores between BZD and OCH A. lyrata tetraploids",
       x='Position along Chromosome 1', y='Fst')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor.x=element_blank(),
        legend.position='none')
##make a threshold
BZD_OCH_Fst$threshold<-BZD_OCH_Fst$`0.193548`>0.55


##load the pure arenosa (KEH) and pure lyrata (LIC) Fst scan
KEH_LIC_Fst<-read_tsv("KEH_LIC_Fst_output.fst")
head(KEH_LIC_Fst)
#plot the Fst scan for these two species
ggplot(KEH_LIC_Fst,aes(x=`1317`,y=`0.000000`,colour=threshold))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=0.5,linetype=2,colour=1)+
  scale_color_brewer(palette='Set2')+
  theme_bw()+
  labs(title="Fst scores between KEH and LIC\nA. Arenosa and A. lyrata tetraploids",
       x='Position along Chromosome 1', y='Fst',colour="Outlier status")+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor.x=element_blank(),
        legend.position='none')

outliers<-KEH_LIC_Fst$`0.000000`>0.50
head(outliers)

KEH_LIC_Fst$threshold<-KEH_LIC_Fst$`0.000000`>0.5
KEH_LIC_Fst

#####################
mean_KEH_LIC_Fst<-mean(KEH_LIC_Fst$`0.000000`)
outliers<-mean_KEH_LIC_Fst+1.5*IQR(KEH_LIC_Fst$`0.000000`)
KEH_LIC_Fst$outliers<-outliers
head(KEH_LIC_Fst)
