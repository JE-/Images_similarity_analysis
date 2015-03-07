##Setup--------------------------------------------------------------------------------
rm(list=ls())
library(neuroim)
library(multicore)
library(FactoMineR)
subjects <- c("S01","S02","S03","S04","S05","S06","S07","S08","S09","S10","S11")

#PARAMETERS:
# Whole-brain ROI?
Whole_brain_ROI = F

#Collapse across runs. But across durations as well?
split_by_duration = F

# pairwiseCorr or pairwiseRV?
pairwise = "pairwiseCorr"
#--------------------------------------------------------------------------------------

pCorr <- ifelse (pairwise =="pairwiseCorr",T,F)
source("C:/Main/OP/7. images study -- analysis using R scripts/analysis/getAnalysisDirectory.R")

# Do this for all subjects
for (subject in subjects){
  print(subject)
  setwd(paste("C:/Main/OP/7. images study -- analysis using R scripts/",subject,sep=""))
  
  pairwiseRV <- function(matlist, center=TRUE, scale=FALSE, ncomp=2) {
    ### creates all pairwise combinations
    pcomb <- t(combn(1:length(matlist), 2))
    
    ## if center and/or scale are TRUE, then center and/or scale columns
    if (center || scale) {
      matlist <- lapply(matlist, scale, center=center, scale=scale)
    }
    
    Rv <- mean(unlist(lapply(1:nrow(pcomb), function(n) {
      i <- pcomb[n,1]
      j <- pcomb[n,2]
      
      M1 <- matlist[[i]]
      M2 <- matlist[[j]]
      
      pca.1 <- prcomp(M1)
      pca.2 <- prcomp(M2)
      
      coeffRV(pca.1$x[,1:ncomp], pca.2$x[,1:ncomp])$rvstd  
    })))
  }
  
  pairwiseCorr <- function(matlist, center=TRUE, scale=FALSE,roi) {
    ### creates all pairwise combinations
    pcomb <- t(combn(1:length(matlist), 2))
    
    ## if center and/or scale are TRUE, then center and/or scale columns
    if (center || scale) {
      matlist <- lapply(matlist, scale, center=center, scale=scale)
    }
    
    corlist <- lapply(1:nrow(pcomb), function(n) {
      #print(n)#je
      i <- pcomb[n,1]
      j <- pcomb[n,2]
      
      M1 <- matlist[[i]]
      M2 <- matlist[[j]]
      cor(t(M1), t(M2))  
    })
    
    #Now you should have a structure corlist which is a list structure of 276 matrices. Each such matrix is 16X16.
    
    ## take average pairwise correlation ('Reduce' is  a fast way to sum elements of a list)
    cor.avg <- Reduce("+", corlist)/length(corlist)
    
    if (Whole_brain_ROI == T){
      write.table(cor.avg, paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",getAnalysisDirectory(Whole_brain_ROI,split_by_duration,pCorr),"/",subject,"_WholeBrainROI_CorrMat.txt",sep=""))
    }else{
      write.table(cor.avg, paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",getAnalysisDirectory(Whole_brain_ROI,split_by_duration,pCorr),"/",subject,"_ROI",roi,"CorrMat.txt",sep=""))
    }
    
    ## compute mean correlation along diagonal
    cor.diag <- mean(diag(cor.avg))
    ## compute mean correlation of off diagonal elements
    cor.off <- mean(c(cor.avg[lower.tri(cor.avg)],
                      cor.avg[upper.tri(cor.avg)]))
    
    cor.diff <- cor.diag - cor.off
    cor.diff
  }

  
  ##-------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## load design table
  ##design <- read.table("rsa_design.txt", header=TRUE)
  # contains inforamtions regarding the bata images of each condition.
  design <- read.csv("rsa_design.txt", header=TRUE, sep=" ") #je
  
  ## load in Harvard-Oxford atlas that has been resampled to match
  ## the image dimensions of our MNI-normalized betas [61, 73, 61]
  harvard <- loadVolume("C:/Main/OP/7. images study -- analysis using R scripts/S01/HarvardOxford_dup_cortsub_thr0_reduced.nii")
  
  ## load in first image and use nonzero values as mask
  mask <- loadVolume("im_nwtstats.nii", 1)
  bvec <- loadVector("im_nwtstats.nii", mask=which(mask != 0))
  ## dim(bvec) is [61  73  61] 384 (because 16 images * 4 durations * 6 sessions = 384)
  
  ## find unique roi numbers (0-112) -- Here '0' is the "background voxel"
  roinums <- sort(unique(as.vector(harvard)))
  
  if (Whole_brain_ROI == T){
    #to take all of the ROIs
    roinums <- -99
  }else {
    ## exclude first element ('0')
    roinums <- roinums[-1]
  }
  
  ## iterate over roi indices in parallel using 'mclapply' (multicore version of standard R function 'lapply')
  ##cor.res <- mclapply(roinums, function(roi) {
  cor.res <- lapply(roinums, function(roi) {
    print(roi)

    if (Whole_brain_ROI == T){
      #to take all voxels, from all ROIs
      roi.idx <- which(harvard != roi & mask != 0)
    }else {
      ## find 1D indices in harvard volume for current roi num
      ## only include indices in harvard mask that intersect mask
      roi.idx <- which(harvard == roi & mask != 0)
    }
    
    ## extract data matrix for indices, will yield an N X P matrix
    ## N = rows = number of betas (i.e. should be 16*4*6 = 384 rows)
    ## P = columns = number of voxels in current roi 
    roimat <- series(bvec, roi.idx) 
    ## type dim(roimat) to get matrix dimensions as a vector
    
    if (split_by_duration == F){
      ## split full matrix by the run number, yielding a list of matrices, 1 per run
      matlist <- split(as.data.frame(roimat), design$run)
    }else if (split_by_duration == T){
      matlist <- split(as.data.frame(roimat),paste(design$duration,design$run))
      # length(matlist) is 24 because 4 (duarions)* 6 (runs) = 24
      # dim(matlist[[1]]) = 16 1925 because 16 (images) and 1925 (i.e. number of voxels in the current roi)
    }
    
    if (pairwise == "pairwiseCorr"){
      ## compute pairwise correlations for each combinations of runs
      cest <- pairwiseCorr(matlist, center=TRUE, scale=FALSE,roi) #je
      #cest is the difference between the mean diagonal to the mean off diagonal elements.
      print(cest)
      #cest is one value -- the difference between the main vs. off diagonal for the average correlations matrix from all 15 combinations of runs (or more comparisons if serated by runs and durations)
      cest
    } else if (pairwise == "pairwiseRV"){
      ## compute pairwise correlations for each combinations of runs
      cest <- pairwiseRV(matlist, center=TRUE, scale=FALSE) #je
      #cest is the difference between the mean diagonal to the mean off diagonal elements.
      print(cest)
      #cest is one value -- the difference between the main vs. off diagonal for the average correlations matrix from all 15 combinations of runs (or more comparisons if serated by runs and durations)
      cest
    }
  })

  if (pairwise == "pairwiseCorr"){
    if (Whole_brain_ROI == T){
      write.table(as.data.frame(unlist(cor.res)),paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",getAnalysisDirectory(Whole_brain_ROI,split_by_duration,pCorr),"/",subject,"_wholeBrainROI_CorrDiff.txt",sep=""))                 
    }else{
      write.table(as.data.frame(unlist(cor.res)),paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",getAnalysisDirectory(Whole_brain_ROI,split_by_duration,pCorr),"/",subject,"_HarvOxROIs_CorrDiff.txt",sep=""))     
      ## create map of correlation values -- 'fill' maps the roinumbers from the atlas to the supplied correlations
      corvol <- fill(harvard, unlist(cor.res))
      writeVolume(corvol,paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",getAnalysisDirectory(Whole_brain_ROI,split_by_duration,pCorr),"/",subject,"_HarvOxROIs_CorrDiff.nii",sep=""))
    }
  }
  if (pairwise == "pairwiseRV"){
    if (Whole_brain_ROI == T){
      write.table(as.data.frame(unlist(cor.res)),paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",getAnalysisDirectory(Whole_brain_ROI,split_by_duration,pCorr),"/",subject,"_wholeBrainROI_RV.txt",sep=""))      
    }else{
      write.table(as.data.frame(unlist(cor.res)),paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",getAnalysisDirectory(Whole_brain_ROI,split_by_duration,pCorr),"/",subject,"_HarvOxROIs_RV.txt",sep=""))     
      ## create map of correlation values -- 'fill' maps the roinumbers from the atlas to the supplied correlations
      corvol <- fill(harvard, unlist(cor.res))
      writeVolume(corvol,paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",getAnalysisDirectory(Whole_brain_ROI,split_by_duration,pCorr),"/",subject,"_HarvOxROIs_RV.nii",sep=""))
    }
  }
}
