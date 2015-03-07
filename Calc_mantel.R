#-----------------------------------------------------------------------
# Compute correlations between brain ROIs (using the "Mantel" function).
# There are 112 ROIs to begin with, but you can decrease the number of ROIs (drop uninformative ROIs), but specifying a threshold.
# For example, a threshold of 0.05 will only include ROIs that code an image significantly across subjects.
#-----------------------------------------------------------------------
rm(list=ls())
setwd("C:/Main/OP/7. images study -- analysis using R scripts/analysis")
source("getROILAbels.R")
library(neuroim)
library(multicore)
library(FactoMineR)
library(ade4)
library(ggplot2)
library(proxy)
library(DistatisR)
library(clickme)
library(MASS)
subjects <- c("S01","S02","S03","S04","S05","S06","S07","S08","S09","S10","S11")

#PARAMETERS:
P_THRESHOLD <- 0.05
ALL_ROIS <- getROILabels()
#------------------------1. Reduce number of ROIs ? ----------------------------------------------------------------------------------------------
x<-as.data.frame(lapply(subjects,function(x){
  read.table(paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/HarvardOxfordROIs_pairwiseCorr_byRunsAndDurations/",x,"_HarvOxROIs_CorrDiff.txt",sep=""),col.names = x)
}))

tTests<-apply(x,1,function(x) {t.test(x,alternative="greater")})

pValues<-0
for (i in 1:112){
  pValues[i]<-tTests[[i]]$p.value
}

numOfROIsSurvivingThreshold <- sum(pValues<P_THRESHOLD)
percentOfROIsTaken <- numOfROIsSurvivingThreshold/112
rois <- ALL_ROIS[pValues<P_THRESHOLD]
#-----------------------2. Compute mantel matrices-----------------------------------------------------------------------------------------------
roi_nums<-1:length(rois)
for (s in 1:length(subjects)){
  
  # load all 112 16X16 matrices for this subject.
  data<-list()
  for (roi in roi_nums){
    data[[roi]]<- as.matrix(read.table(paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",subjects[s],"_ROI",roi,"CorrMat.txt",sep="")))#,col.names=1:16,),rownames.force=FALSE)      
  }
  
  # convert correlation matrices into distance matrices
  for (roi in 1:length(rois)){
    mat<-data[[roi]] 
    newMat<-matrix(NA,16,16)
    for (i in 1:16){
      for (j in 1:16){
        if (i<j){
          newMat[j,i]<-(mat[i,j]+mat[j,i])/2
        }
      }
    }
    newMat <- 1-newMat
    data[[roi]]<-newMat
  }
  
  obsMat<-matrix(NA,length(rois),length(rois))
  mantel_pValues<-matrix(NA,length(rois),length(rois))
  for (i in 1:length(rois)){
    for (j in 1:length(rois)){
      observation <- mantel.rtest(as.dist(data[[i]]), as.dist(data[[j]])) 
      print(i)
      print(observation$obs)
      obsMat[i,j]<-observation$obs
      mantel_pValues[i,j] <- observation$pvalue
    }
  }
  write.table(obsMat, paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",subjects[s],"_obsMat_ReducedROIs.txt",sep=""))
  write.table(mantel_pValues, paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",subjects[s],"_mantel_pValues_ReducedROIs.txt",sep=""))
}
#------------------------------ Read Mantel Matrices -------------------------------------------------------------------------------------------------

subjects <- c("S01","S02","S03","S04","S05","S06","S07","S08","S09","S10","S11")

x<-(lapply(subjects,function(x){
  read.table(paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",x,"_obsMat_ReducedROIs.txt",sep=""))
}))

# Calculate an average observation matrix
obsMatSum<-matrix(0,length(rois),length(rois))
for (i in 1:length(subjects)) {
  obsMatSum <-  obsMatSum + x[[i]]
}
obsMatSum<-obsMatSum/length(subjects)

#--------------------------- plot-------------------------------------------------------------------------------------
new.palette <- colorRampPalette(c("black","red","yellow","white"),space="rgb") 
m <- as.matrix(obsMatSum)
m <- m[length(rois):1,]
# colnames(m) <- 1:112
# rownames(m) <-112:1
colnames(m) <- rois
rownames(m) <-rev(rois)
levelplot(m,scales=list(x=list(rot=90)),main="mean Mantel R correlations between ROIs", xlab="", ylab="", col.regions=new.palette(200), at=seq(min(m,na.rm="T"),max(m,na.rm="T"),0.01))
#-------------------------------------------------------------------------------------------------------------

#-----------------------------------DISTATIS--------------------------------------------------------------------------
#load("threshold.Rdata")
a <- array(dim=c(length(rois),length(rois),11))
for (i in 1:11){
  x <- read.table(paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",subjects[i],"_obsMat_ReducedROIs.txt",sep=""))
  a[,,i] <- as.matrix(x)
}

a <- 1-a #convert to distance
rownames(a) <- rois
distatis_a <- distatis(a)

# 5.2 a compromise plot
factor1 <- 1
factor2 <- 2

#labels <- getROILabels()
compromise.graph.out <- GraphDistatisCompromise(distatis_a$res4Splus$F,axis1 = factor1, axis2 = factor2, constraints = NULL, item.colors = NULL, 
                                                ZeTitle = "Distatis-Compromise (Whole Brain)", nude = FALSE, Ctr = NULL)




# #plot
# new.palette <- colorRampPalette(c("black","red","yellow","white"),space="rgb") 
# colnames(tTests.tValue) <- c(1:112)
# levelplot(tTests.tValue,scales=list(x=list(rot=90)),main="mean Mantel R correlations between ROIs", xlab="ROI Number", ylab="ROI Number", col.regions=new.palette(5000), at=seq(min(tTests.tValue,na.rm="T"),max(tTests.tValue,na.rm="T"),0.1))
# #-------------------------------------------------------------------------------------------------------------


#---------------------------------isoMDS--------------------->
library(MASS)
d <- dist(obsMatSum) # euclidean distances between the rows
mds <- isoMDS(d, k=2) # k is the number of dim
mds # view results

# plot solution 
x <- mds$points[,1]
y <- mds$points[,2]
plot(x, y, xlab="1", ylab="2", 
     main="MDS", type="n")
text(x, y, labels = rois, cex=.8)
#------------------------------------------------------------>






# ------- CODE NOT USED------------------>
# --------------------------- PCA on the average mantel matrix --------------------- ?
# m <- as.matrix(obsMatSum)
# diag(m) <- 1
# for (i in 1:nrow(m)){
#   for (j in 1:ncol(m)){
#     if (i>j)
#       m[i,j]<-m[j,i]
#   }
# }
# pca.m <- prcomp(m,center=T,scale=T)
# summary(pca.m)
# screeplot(pca.m,main="Scree plot",xlab="Components")
#-------------------------------------------------------------------------------------

# # T-test on mantel values to see if greater than 0.
# tTests.tValue<-matrix(NA,112,112)
# tTests.pValue<-matrix(NA,112,112)
# for (i in 1:111){
#   for (j in (i+1):112){
#     print(i)
#     a <- rep(NA, 11)
#     for (s in 1:11){
#       a[s]<-x[[s]][i,j]
#     }
#     ttest <- t.test(a,alternative="greater")
#     tTests.tValue[i,j] <- ttest$statistic
#     tTests.pValue[i,j] <- ttest$p.value
#   }
# }
# 
# #What percentage of ROIs did we take?
# sum(tTests.pValue<P_THRESHOLD,na.rm=T)/(112*112-112)
