getAnalysisDirectory <- function(Whole_brain_ROI, split_by_duration, pairwiseCorr) {
  if (Whole_brain_ROI == T && split_by_duration == T && pairwiseCorr == T)
    "WholeBrain_pairwiseCorr_byRunsAndDurations"
  else if (Whole_brain_ROI == T && split_by_duration == T && pairwiseCorr == F)
    "WholeBrain_pairwiseRV_byRunsAndDurations"
  else if (Whole_brain_ROI == T && split_by_duration == F && pairwiseCorr == T)
    "WholeBrain_pairwiseCorr_byRuns" 
  else if (Whole_brain_ROI == T && split_by_duration == F && pairwiseCorr == F)
    "WholeBrain_pairwiseRV_byRuns"
  else if (Whole_brain_ROI == F && split_by_duration == T && pairwiseCorr == T)
    "HarvardOxfordROIs_pairwiseCorr_byRunsAndDurations"
  else if (Whole_brain_ROI == F && split_by_duration == T && pairwiseCorr == F)
    "HarvardOxfordROIs_pairwiseRV_byRunsAndDurations"
  else if (Whole_brain_ROI == F && split_by_duration == F && pairwiseCorr == T)
    "HarvardOxfordROIs_pairwiseCorr_byRuns"
  else if (Whole_brain_ROI == F && split_by_duration == F && pairwiseCorr == F)
    "HarvardOxfordROIs_pairwiseRV_byRuns"
}