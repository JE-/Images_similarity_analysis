#-----------------------
# Take the 16x16 matrix of each subject (11 such matrices in total) and compute distatis
#-----------------------
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
library(gridExtra)
library(jpeg)

subjects <- c("S01","S02","S03","S04","S05","S06","S07","S08","S09","S10","S11")

a <- array(dim=c(16,16,11))

for (i in 1:11){
  x <- read.table(paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/",subjects[i],"_WholeBrainROI_CorrMat.txt",sep=""))
  a[,,i] <- as.matrix(x)
}

# distance measure instead of a correlation
a <- 1-a
rownames(a) <- c("Beach","Candy","Crying_Baby","Delacroix_Painting","Eyes","Fire","Food","Jordan","Kandinsky_Painting","Massage","Shark","Skyscrapers","Thorns","Tools","Torso","Zebras")

#read images
image_names <-  c("Beach","Candy","Crying_Baby","Delacroix_Painting","Eyes","Fire","Food","Jordan","Kandinsky_Painting","Massage","Shark","Skyscrapers","Thorns","Tools","Torso","Zebras")
image_names <- paste(image_names,".JPG",sep="")

distatis_a <- distatis(a)
# 4.1 Get the bootstrap factor scores (with default 1000 iterations)
#BootF <- BootFactorScores(distatis_a$res4Splus$PartialF)

# 4.2 Get the boostrap from full bootstrap (default niter = 1000)
#F_fullBoot <- BootFromCompromise(a,niter=1000)

#-----------------------------------------------------------------------------
# 5. Create the Graphics
# An RV map
#rv.graph.out <- GraphDistatisRv(distatis_b$res4Cmat$G)

# factors <- as.data.frame(distatis_b$res4Splus$F)
# colnames(factors)<-c("Factor_1","Factor_2","Facotr_3")
# g <- ggplot(factors,aes(Factor_1,Factor_2))
# p <- g + geom_point()
# print(p)

# plot factors 1 and 2
factor1 <- 2
factor2 <- 3

factors <- as.data.frame(distatis_a$res4Splus$F)

xMin <- min(factors[,factor1])
xMax <- max(factors[,factor1])
yMin <- min(factors[,factor2])
yMax <- max(factors[,factor2])

g <- list(dim=c(16))
for (i in 1:16){
  print(i)
  img <- readJPEG(image_names[i])
  g[[i]] <- rasterGrob(img, interpolate=TRUE)
}

scale <- 0.08
qplot(factors[,factor1],factors[,factor2], geom="blank") +
#   for (i in 1:16){
#     annotation_custom(g[[i]], xmin=factors[i,1], xmax=factors[i,1]+scale, ymin=factors[i,2], ymax=factors[i,2]+scale)+
#   }
  annotation_custom(g[[1]], xmin=factors[1,factor1]-scale/2, xmax=factors[1,factor1]+scale/2, ymin=factors[1,factor2]-scale/2, ymax=factors[1,factor2]+scale/2)+
  annotation_custom(g[[2]], xmin=factors[2,factor1]-scale/2, xmax=factors[2,factor1]+scale/2, ymin=factors[2,factor2]-scale/2, ymax=factors[2,factor2]+scale/2)+
  annotation_custom(g[[3]], xmin=factors[3,factor1]-scale/2, xmax=factors[3,factor1]+scale/2, ymin=factors[3,factor2]-scale/2, ymax=factors[3,factor2]+scale/2)+
  annotation_custom(g[[4]], xmin=factors[4,factor1]-scale/2, xmax=factors[4,factor1]+scale/2, ymin=factors[4,factor2]-scale/2, ymax=factors[4,factor2]+scale/2)+
  annotation_custom(g[[5]], xmin=factors[5,factor1]-scale/2, xmax=factors[5,factor1]+scale/2, ymin=factors[5,factor2]-scale/2, ymax=factors[5,factor2]+scale/2)+
  annotation_custom(g[[6]], xmin=factors[6,factor1]-scale/2, xmax=factors[6,factor1]+scale/2, ymin=factors[6,factor2]-scale/2, ymax=factors[6,factor2]+scale/2)+
  annotation_custom(g[[7]], xmin=factors[7,factor1]-scale/2, xmax=factors[7,factor1]+scale/2, ymin=factors[7,factor2]-scale/2, ymax=factors[7,factor2]+scale/2)+
  annotation_custom(g[[8]], xmin=factors[8,factor1]-scale/2, xmax=factors[8,factor1]+scale/2, ymin=factors[8,factor2]-scale/2, ymax=factors[8,factor2]+scale/2)+
  annotation_custom(g[[9]], xmin=factors[9,factor1]-scale/2, xmax=factors[9,factor1]+scale/2, ymin=factors[9,factor2]-scale/2, ymax=factors[9,factor2]+scale/2)+
  annotation_custom(g[[10]],xmin=factors[10,factor1]-scale/2,xmax=factors[10,factor1]+scale/2,ymin=factors[10,factor2]-scale/2, ymax=factors[10,factor2]+scale/2)+
  annotation_custom(g[[11]],xmin=factors[11,factor1]-scale/2,xmax=factors[11,factor1]+scale/2,ymin=factors[11,factor2]-scale/2, ymax=factors[11,factor2]+scale/2)+
  annotation_custom(g[[12]],xmin=factors[12,factor1]-scale/2,xmax=factors[12,factor1]+scale/2,ymin=factors[12,factor2]-scale/2, ymax=factors[12,factor2]+scale/2)+
  annotation_custom(g[[13]],xmin=factors[13,factor1]-scale/2,xmax=factors[13,factor1]+scale/2,ymin=factors[13,factor2]-scale/2, ymax=factors[13,factor2]+scale/2)+
  annotation_custom(g[[14]],xmin=factors[14,factor1]-scale/2,xmax=factors[14,factor1]+scale/2,ymin=factors[14,factor2]-scale/2, ymax=factors[14,factor2]+scale/2)+
  annotation_custom(g[[15]],xmin=factors[15,factor1]-scale/2,xmax=factors[15,factor1]+scale/2,ymin=factors[15,factor2]-scale/2, ymax=factors[15,factor2]+scale/2)+
  annotation_custom(g[[16]],xmin=factors[16,factor1]-scale/2,xmax=factors[16,factor1]+scale/2,ymin=factors[16,factor2]-scale/2, ymax=factors[16,factor2]+scale/2)+
  labs(title = paste("Distatis: Factors ",factor1," and ",factor2,sep=""))
  
# 5.2 a compromise plot
compromise.graph.out <- GraphDistatisCompromise(distatis_a$res4Splus$F,axis1 = factor1, axis2 = factor2, constraints = NULL, item.colors = NULL, 
                                                                                 ZeTitle = "Distatis-Compromise (Whole Brain)", nude = FALSE, Ctr = NULL)


# 5.3 a partial factor score plot
partial.scores.graph.out <-GraphDistatisPartial(distatis_a$res4Splus$F,distatis_a$res4Splus$PartialF)

# 5.4 a bootstrap confidence interval plot
#5.4.1 with ellipses
boot.graph.out.ell <- GraphDistatisBoot(distatis_a$res4Splus$F,BootF)
#or
# boot.graph.out <- GraphDistatisBoot(testDistatis$res4Splus$F,F_fullBoot)
#5.4.2 with hulls
boot.graph.out.hull <- GraphDistatisBoot(testDistatis$res4Splus$F,BootF,ellipses=FALSE)

