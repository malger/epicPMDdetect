#' Estimates AlphaValues based on window offsets for a Genomic Range
#' @param m Granges including Total,Methylated Counts
#' @param offsets list of offsets describing the datapoints to be used to calc the alpha values for each parameter
#' @return estimated alpha values
#' @description The list offsets musst be in the format:
#'  list(offsetConfig1->list(chr1,2,..),offsetConfig2->list(chr1,2,..))
#'  Offsets are relative to the datapoint position. (-1,-2,-3,0,1,2,3)
#' @export
alphaValueEstimation <-
  function(m,offsets){
    
    chrvec = seqnames(m)
    chrs = names(offsets[[1]])
    meth = m$M
    total = m$`T`
    
    alphas = foreach(koffsets = offsets,.packages = c('S4Vectors','foreach','zoo','GenomicRanges'),.export=c("as.boolean","calculateAlphaDistr"),.noexport = c('m','offsets')) %:% 
      foreach(koff.chr = koffsets,chr.sel=chrs,.final = function(x)setNames(x,chrs)) %dopar% {
        print(chr.sel)
        chrindx=as.boolean(chrvec==chr.sel)
        calculateAlphaDistr(meth[chrindx],  total[chrindx] ,koff.chr)
      }
    names(alphas)  = names(offsets)
    alphas
    
}
