#' Generates a genome segmentation including PMDs,nonPMDs and their unsure counterparts.'
#' @param m GRange representing the EPIC datapoints. Metadata musst include Total,Methylated Counts
#' @param y.list prediction object of library(mhsmm).
#' @param seqLengths Seqlenghts objects representing the assembly sequence lengths of the chromosomes.
#' @return GRange genome segmentation with PMDs,nonPMDs and their unsure counterparts.
createGRangesObjectPMDSegmentation <-
function(m.list, y.list, seqLengths){

  message("creating GRanges object")

 
        
  segList <- foreach(chr.sel = unique(seqnames(m)),
                     m.chr= m.list,
                     y = y.list,
                     .noexport = c("m","y.list"),
                     .export = c("seqLengths"),
                     .packages=c('GenomicRanges','IRanges')) %dopar% {

    n = length(m.chr)
    segCG = y$s
                       
    mids <- round(0.5*(start(m.chr)[-length(m.chr)] + start(m.chr)[-1])) #get mids between dps
   
    #copy dp assignment (PMD/notPMD/unkown) to Rle having the bp distance as runLength
    segNt <- Rle(lengths=c(diff(c(1,mids)),seqLengths[chr.sel]-mids[n-1]+1), values=as.numeric(segCG)) 
    #convert to Genomic Ranges
    segChr <- GRanges(seqnames=chr.sel, 
                      IRanges(start=c(1,cumsum(runLength(segNt))[-nrun(segNt)]+1), 
                              end=cumsum(runLength(segNt))), 
                      strand="*", 
                      type=c("notPMD","PMD","unkown")[runValue(segNt)], 
                      nCG=runLength(segCG), 
                      seqlengths=seqLengths
              )
    segChr

  }

  segments <- do.call(c, unname(segList))
  segments

}
