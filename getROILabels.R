getROILabels <- function() {
  HO <- read.table(paste("C:/Main/OP/7. images study -- analysis using R scripts/analysis/Harvard_labels.txt"),header=T,sep="")
  HO <- HO[order(HO$ROI),]
  rois <- array()
  
  for (i in 1:112){
    if (i %in% HO$ROI){           
      rois[i]<-paste(ifelse(HO[which(HO$ROI==i),5]=="left","L.","R."),HO[which(HO$ROI==i),2])
    } else{
      rois[i]<-"NA"
    }
  } 
  rois
}
