dynSegmentPMDs <-
function(m,
         chr.sel,
         pdfFilename = NULL,
         seqLengths = "hg19",
         num.cores = 4,
         method = 'knn',
         tilewidth = 2.5*10^5,
         k=30,
         q=1,
         cutoff = NULL
        ){ 

        bins = NULL
        if(seqLengths == "hg19"){
                library(TxDb.Hsapiens.UCSC.hg19.knownGene)
                seqLengths = seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)
        }
        seqLengths = filterPresentChromosomesNSort(seqLengths,m)
        if(method == 'tiles'){
          bins <- getDynamicCpGBins(m,seqLengths,tileWidth)
        }
        hmm.model <- trainPMDHMM(m, chr.sel, num.cores, plot.distr=TRUE, pdfFilename,method,bins,k,q,cutoff)
        y.list <- PMDviterbiSegmentation(m, hmm.model, num.cores,method,bins,k,q,cutoff)
        segments <- createGRangesObjectPMDSegmentation(m, y.list, num.cores, seqLengths)
        segments
         
 }


filterPresentChromosomesNSort = function(seqlengths,gr){
        seqlengths = seqlengths[names(seqlengths) %in% levels(seqnames(gr))]
        
        seqlengths[order(match(names(seqlengths),names(seqlengths(gr))))]
}
