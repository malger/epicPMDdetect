dynSegmentPMDs <-
function(m,
         chr.sel,
         pdfFilename = NULL,
         seqLengths,
         num.cores = 4,
         method = 'knn',
         tilewidth = 2.5*10^5,
         k=30,
         q=1
        ){ 

        bins = NULL
        if(method == 'tiles'){
          bins <- getDynamicCpGBins(m,tileWidth)
        }
        hmm.model <- trainPMDHMM(m, chr.sel, num.cores, plot.distr=TRUE, pdfFilename,method,bins,k,q)
        y.list <- PMDviterbiSegmentation(m, hmm.model, num.cores,method,bins,k,q)
        segments <- createGRangesObjectPMDSegmentation(m, y.list, num.cores, seqlengths(m))
        segments
         
 }
