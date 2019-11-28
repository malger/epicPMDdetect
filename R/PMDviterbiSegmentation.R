PMDviterbiSegmentation <-
function(m, hmm.model, num.cores,method,bins,k,q,cutoff = NULL){

  low_density_cutoff = 0.2
  message("performing viterbi segmentation")


  y.list=mclapply(unique(seqnames(m)), function(chr.sel,mc.cores=num.cores){
  
    chrindx=as.character(GenomicRanges::seqnames(m))==chr.sel;
    Total <- as.numeric(values(m[chrindx])[, 1])
    Meth <- as.numeric(values(m[chrindx])[, 2])
    
    dist = start(m[chrindx][-1])-end(m[chrindx][-length(m[chrindx])])
    
    offsets = 101
    if(method == 'knn')
      offsets = getKNNOffsets(GenomicRanges::start(m[chrindx]),k,q,cutoff)
    if (method == 'tiles')
      offsets = getBinSizeOffsets(bins[[chr.sel]],q)
    
    score <- calculateAlphaDistr(Meth, Total, num.cores,offsets)
    train=list(x=score, N=length(score));
    y=predict(hmm.model, train);

    #remove regions that are too short
    ttt=Rle(y$s)
    
    if(length(runValue(ttt))>1){
      # first mark PMDs, that are based on low density points
      pmdsIndx= runValue(ttt)==2 #any PMD 
      pmdsStart = cumsum(runLength(ttt))[pmdsIndx]
      pmdProp= lapply(1:length(pmdsStart),function(i){
        range = na.omit(dist[pmdsStart[i]:(pmdsStart[i]+runLength(ttt)[pmdsIndx][i])])
        c(length = sum(range),Ndps = length(range),density = length(range)/sum(range))
      })
      df <- data.frame(matrix(unlist(pmdProp), nrow=length(pmdProp), byrow=T))
      colnames(df) = names(pmdProp[[1]])
      lq = quantile(df$density,probs=low_density_cutoff,na.rm=TRUE)
      runValue(ttt)[pmdsIndx][df$density<lq]=3 #low PMD point coverage
      
      NpmdsIndx = runValue(ttt) == 1
      NpmdsStart = cumsum(runLength(ttt))[NpmdsIndx]
      NpmdProp= lapply(1:length(NpmdsStart),function(i){
        range = na.omit(dist[NpmdsStart[i]:(NpmdsStart[i]+runLength(ttt)[NpmdsIndx][i])])
        c(length = sum(range),Ndps = length(range),density = length(range)/sum(range))
      })
      df <- data.frame(matrix(unlist(NpmdProp), nrow=length(NpmdProp), byrow=T))
      print(head(df))
      colnames(df) = names(NpmdProp[[1]])
      lq = quantile(df$density,probs=low_density_cutoff,na.rm=TRUE)
      runValue(ttt)[NpmdsIndx][df$density<lq]=4 #low nPMD point coverage
      
    }
      
    
      
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
    
    y$s=as.vector(ttt);
    y$score=score
    print(paste0('finished chromosome:',chr.sel))
    
    return(y)

  }, mc.cores=num.cores);
  names(y.list) = unique(seqnames(m))
  y.list

}
