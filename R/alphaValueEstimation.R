alphaValueEstimation <-
  function(m,offsets){
    
    message("estimating alpha values with k=",k)
    
    
    y.list=foreach(chr.sel = unique(seqnames(m))) %dopar% {
      
      chrindx=as.character(GenomicRanges::seqnames(m))==chr.sel;
      Total <- as.numeric(values(m[chrindx])[, 1])
      Meth <- as.numeric(values(m[chrindx])[, 2])
      score <- calculateAlphaDistr(Meth, Total ,offsets[[chr.sel]])

      return(score)
      
    }
    names(y.list) = unique(seqnames(m))
    y.list
    
  }
