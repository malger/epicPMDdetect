#' Predicts PMDs/nPMD segments based on HiddenMarkovModel. Selection by ML using viterbi selection.
#' @param alphas estimated alpha values for each datapoint. format:list(chr->vector(alphaValues))
#' @param hmm.model trained Hidden Markov Model
#' @param markHighDPdistVarSegms mark segments with high variance in pairwise dp distance
#' @param knns list(k->list(chrs->FNN::get.knnx objs)), only needed if markHighDPdistVarSegms is TRUE
#' @return prediction of PMD/nPMD segments for each chromosome
#' @export
PMDviterbiSegmentation <-
function(alphas, hmm.model,markHighDPdistVarSegms=F,knns=NULL){

  message("performing viterbi segmentation")


  y.list = foreach(chr.sel = names(alphas),.export = 'getHighVariableDistNNs') %dopar%{
  
    message('predicting chromosome ',chr.sel)
    score = alphas[[chr.sel]]
    train=list(x=score, N=length(score));
    y=predict(hmm.model, train);
    
    ttta = rle(y$s)

    if(markHighDPdistVarSegms == T){
      
      outls = which(getHighVariableDistNNs(knns[[chr.sel]]))
      tt = with(ttta, data.frame(number = values,
                                     start = cumsum(lengths) - lengths + 1,
                                     end = cumsum(lengths))[order(values),])
      
      if(length(outls)>0){
        i = outls[1]
        while(i<=length(y$s)){
          j = which(i>=tt$start & i<=tt$end)
          print(j)
          
          ttta$values[j] = ifelse( ttta$values[j] == 2,3,4)
          blockend = tt[j,"end"]
          i = outls[min(which(outls>blockend))]
          print(i)
          if(is.na(i)){
            break;
          }
        }
      }
    }
      
      # if(length(runValue(ttt))>1){
      #   # first mark PMDs, that are based on low density points
      #   pmdsIndx= runValue(ttt)==2 #any PMD 
      #   pmdsStart = cumsum(runLength(ttt))[pmdsIndx]
      #   pmdProp= lapply(1:length(pmdsStart),function(i){
      #     range = na.omit(dist[pmdsStart[i]:(pmdsStart[i]+runLength(ttt)[pmdsIndx][i])])
      #     c(length = sum(range),Ndps = length(range),density = length(range)/sum(range))
      #   })
      #   df <- data.frame(matrix(unlist(pmdProp), nrow=length(pmdProp), byrow=T))
      #   colnames(df) = names(pmdProp[[1]])
      #   lq = quantile(df$density,probs=low_density_cutoff,na.rm=TRUE)
      #   runValue(ttt)[pmdsIndx][df$density<lq]=3 #low PMD point coverage
      #   
      #   NpmdsIndx = runValue(ttt) == 1
      #   NpmdsStart = cumsum(runLength(ttt))[NpmdsIndx]
      #   NpmdProp= lapply(1:length(NpmdsStart),function(i){
      #     range = na.omit(dist[NpmdsStart[i]:(NpmdsStart[i]+runLength(ttt)[NpmdsIndx][i])])
      #     c(length = sum(range),Ndps = length(range),density = length(range)/sum(range))
      #   })
      #   df <- data.frame(matrix(unlist(NpmdProp), nrow=length(NpmdProp), byrow=T))
      #   print(head(df))
      #   colnames(df) = names(NpmdProp[[1]])
      #   lq = quantile(df$density,probs=low_density_cutoff,na.rm=TRUE)
      #   runValue(ttt)[NpmdsIndx][df$density<lq]=4 #low nPMD point coverage
      #   
      # }
      # 
   

    
      
    # # first take regions that are PMD, but too short and make them nonPMD
    # min.length = round(quantile(runLength(ttt)[runValue(ttt==2)]))[3]
    # indx=runLength(ttt)<min.length & runValue(ttt)==2 #any PMD that is shorter 2
    # pos = cumsum(runLength(ttt))[indx]
    # remL = sapply(1:length(pos),function(i){
    #   range = na.omit(dist[pos[i]:(pos[i]+runLength(ttt)[indx][i])])
    #   density = (length(range)/sum(range))
    #   
    #   if (sum(range) < 50000){
    #     1
    #   }else{
    #     
    #   
    #     if(density>0.0005){
    #       0
    #     } else {
    #       1
    #     }
    #   }
    # })
    # runLength(ttt)[indx][remL] = 1 
    
    # now vice versa
    
    # min.length = round(quantile(runLength(ttt)[runValue(ttt==1)]))[3]
    # indx=runLength(ttt)<=min.length & runValue(ttt)==1;
    # pos = cumsum(runLength(ttt))[indx]
    # 
    # remL = sapply(1:length(pos),function(i){
    #   print(i)
    #   range = na.omit(dist[pos[i]:(pos[i]+runLength(ttt)[indx][i])])
    #   density = (length(range)/sum(range))
    #   # print(paste(sum(range),density))
    #     if(density>0.0005){
    #       0
    #     } else {
    #       1
    #     }
    #      
    #   
    #   
    # })
    # runValue(ttt)[indx2][remL]=2;
    
    y$s=inverse.rle(ttta)
    y$score=score
    print(paste0('finished chromosome:',chr.sel))
    
    return(y)

  }
  names(y.list) = names(alphas)
  y.list

}
