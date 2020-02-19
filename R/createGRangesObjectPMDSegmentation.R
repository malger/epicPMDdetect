#' Generates a genome segmentation including PMDs,nonPMDs and their unsure counterparts.'
#' @param m GRange representing the EPIC datapoints. Metadata musst include Total,Methylated Counts
#' @param y.list prediction object of library(mhsmm).
#' @param seqLengths Seqlenghts objects representing the assembly sequence lengths of the chromosomes.
#' @return GRange genome segmentation with PMDs,nonPMDs and their unsure counterparts.
createGRangesObjectPMDSegmentation <-
function(m, y.list, seqLengths){

  message("creating GRanges object")
  chrs=names(y.list)

  segList <- foreach(chr.sel = chrs,.packages=c('GenomicRanges','IRanges')) %dopar% {

    indx=as.character(seqnames(m))==chr.sel
    n <- sum(indx)
    mids <- round(0.5*(start(m[indx])[-length(m[indx])] + start(m[indx])[-1]))
    segCG <- Rle(y.list[[chr.sel]]$s)
    artifical.ylist.ext = c(y.list[[chr.sel]]$s,
                            rep(tail(y.list[[chr.sel]]$s,1),length(m[indx])-length(y.list[[chr.sel]]$s))
    )
    
    segNt <- Rle(lengths=c(diff(c(1,mids)),seqLengths[chr.sel]-mids[n-1]+1), values=artifical.ylist.ext)
    segChr <- GRanges(seqnames=chr.sel, 
                      IRanges(start=c(1,cumsum(runLength(segNt))[-nrun(segNt)]+1), 
                              end=cumsum(runLength(segNt))), 
                      strand="*", 
                      type=c("notPMD","PMD","likelyPMD","likelynotPMD")[runValue(segNt)], 
                      nCG=runLength(segCG), 
                      seqlengths=seqLengths
              )
    segChr

  }

  segments <- do.call(c, unname(segList))
  segments

}
