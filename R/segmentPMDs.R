segmentPMDsKNN <-
function(m,
         chr.sel="chr2",
         pdfFilename = NULL,
         seqLengths = "hg19",
         method = 'knn',
         knn=c(20,25,30),
         q=1,
         cutoff = 50000,
         num.cores = 1,
         markLowDensitySegm = T
        ){ 

        parallel.setup(num.cores)
        if(seqLengths == "hg19"){
                seqLengths = seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)
        }
        seqLengths = filterPresentChromosomesNSort(seqLengths,m)
        
        offsets = calculateKNNs(m,knn,cutoff,q)
        alphas = foreach(k = names(offsets)) %dopar%  alphaValueEstimation(m,offsets[[k]])
        
        names(alphas) = paste0('knn',knn)
        
        df = do.call(cbind,alphas)
        comb_alpha = foreach(i = 1:nrow(df)) %dopar%{
                apply(sapply(df[i,],cbind),1,mean) #TODO remove hack
        }
        #comb_alpha = lapply(comb_alpha,'/',length(knn))
        names(comb_alpha) = names(alphas[[1]])
        
        
        hmm.model <- trainPMDHMM(comb_alpha, chr.sel, num.cores, plot.distr=TRUE, pdfFilename,k,q,cutoff)
        y.list <- PMDviterbiSegmentation(comb_alpha, hmm.model, num.cores,k,q,cutoff,markLowDensitySegm)
        segments <- createGRangesObjectPMDSegmentation(m, y.list, num.cores, seqLengths)
        segments
         
 }


filterPresentChromosomesNSort = function(seqlengths,gr){
        seqlengths = seqlengths[names(seqlengths) %in% levels(seqnames(gr))]
        seqlengths[order(match(names(seqlengths),names(seqlengths(gr))))]
}
