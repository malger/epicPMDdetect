#' Predicts PMDs/nPMD segments based on HiddenMarkovModel. Selection by ML using viterbi selection.
#' @param alphas estimated alpha values for each datapoint. format:list(chr->vector(alphaValues))
#' @param hmm.model trained Hidden Markov Model
#' @param markHighDPdistVarSegms mark segments with high variance in pairwise dp distance
#' @param knns list(k->list(chrs->FNN::get.knnx objs)), only needed if markHighDPdistVarSegms is TRUE
#' @return prediction of PMD/nPMD segments for each chromosome
#' @export
PMDviterbiSegmentation <-
function(alphas, hmm.model,markLowDensitySegm=T,m.list){

  message("performing viterbi segmentation")
  
  #calc next neighbor distance
  ndp_dist = start(m)[-1]-start(m)[-length(m)]
  ndp_dist[ndp_dist<0] = NA # change of chromosome fix
  m$nndist = c(ndp_dist,NA)
  m$critical = rep(F)
 
  

  y.list = foreach(chr.sel = names(alphas),.export=c("Rle")) %dopar%{
  
    message('predicting chromosome ',chr.sel)
    score = alphas[[chr.sel]]
    train=list(x=score, N=length(score));
    y=predict(hmm.model, train);
    
    segms = Rle(y$s)

    if (length(segms)>1 && markLowDensitySegm){
      
      #calc distribution quantiles , one is interested in dps with exceptional high distance to the next
      high99 = quantile(m[seqnames(m)==chr.sel]$nndist,na.rm=T,probs =0.995)
      # mark data points that are in the >99% area
      cdpsPos = na.omit(which(m[seqnames(m)==chr.sel]$nndist>high99))
      max_val = tail(which(seqnames(m)==chr.sel),1)
      cdpsPos = c(cdpsPos,sapply(cdpsPos+1,min,max_val))
      segms[cdpsPos] = 3
      
   
      #  mark PMDs, that are based on data points with high next DP distance
      # #chromosome criticals
      # cdpsPos= which(m[seqnames(m)==chr.sel]$critical)
      # 
      # #selector for rle elements that are pmd
      # pmdIndx = as.logical(runValue(segms)==2)
      # #get the end positions of the pmd segments
      # endpos = cumsum(runLength(segms))
      # #convert to start points
      # startPos = c(1,head(endpos+1,-1))#
      # #filter pmds
      # pmdsStartPos= startPos[pmdIndx] 
      # 
      # 
      # #find nearest start position of the criticals, to determin their parent segments
      # 
      # # this functions finds the indx of the next pos. integer that minimizes the distance
      # # thereby we get the start positions
      # nearest = function(num,arr) {
      #   diffs = num-arr
      #   diffs = diffs[diffs>0]
      #   which.min(diffs)
      # }
      # # apply for all critical dps, get their parent PMD blocks
      # criticalPMDsIndx = unlist(lapply(cdpsPos,nearest,pmdsStartPos))
      # #get all unique pmds
      # affected_pmds = unique(criticalPMDsIndx)
      # #get the number of dps inside pmd 
      # pmdsNdps= (runLength(segms))[pmdIndx][affected_pmds]
      # #count how many points inside pmd are critical
      # nCriticalDPsPerPMD = table(criticalPMDsIndx)
      # #get the propotion of critical dps for each affected pmd
      # propCritical_nDPs = nCriticalDPsPerPMD/pmdsNdps
      # 
      # markPMDs = affected_pmds[which(propCritical_nDPs>.1)] #mark pmds that have more than 10% criticals
      # runValue(segms)[pmdIndx][affected_pmds] = 3
      
    }
        
    
    y$s=segms
    y$score=score
    #print(paste0('finished chromosome:',chr.sel))
    
    return(y)

  }
  names(y.list) = names(alphas)
  y.list

}
