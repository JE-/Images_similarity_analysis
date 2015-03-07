# compute a one-tailed t test for each ROI (with all d values from all 11 participants) to test if that ROI successfully coded the same image vs. different images.

rm(list=ls())
setwd("C:/Main/OP/7. images study -- analysis using R scripts/analysis")
library(neuroim)
library(multicore)
library(FactoMineR)
library(compute.es)
library(ade4)
library(ggplot2)
library(proxy)
library(DistatisR)

subjects <- c("S01","S02","S03","S04","S05","S06","S07","S08","S09","S10","S11")

x<-as.data.frame(lapply(subjects,function(x){
  read.table(paste("C:/Main/OP/7. images study -- analysis using R scripts/",x,"/harvard_cor_diag_byRunsAndDurations.txt",sep=""),col.names = x)
}))

tTests<-apply(x,1,function(x) {t.test(x,alternative="greater")})

# pValues<-0
# for (i in 1:112)
#   if (tTests[[i]]$p.value < 0.05)
#     pValues[i]<-1
#   else
#     pValues[i]<-0
# }

harvard <- loadVolume("C:/Main/OP/7. images study -- analysis using R scripts/S01/HarvardOxford_dup_cortsub_thr0_reduced.nii")

pValues<-0
for (i in 1:112){
  pValues[i]<-tTests[[i]]$p.value
}

tValues<-0
for (i in 1:112){
  tValues[i]<-tTests[[i]]$statistic
}

pValuesBrain <- fill(harvard, pValues)
writeVolume(pValuesBrain, "pValuesBrain_byRunsAndDurations.nii")

tValuesBrain <- fill(harvard, tValues)
writeVolume(tValuesBrain, "tValuesBrain_byRunsAndDurations.nii")
              
# z<-data.frame(rowMeans(x)>0,pValues<0.05)
