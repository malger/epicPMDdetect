#' Segments a genome into PMD and nonPMD regions using an HMM and KNNs
#' 
#' @param m Granges including Metadata for Total,Methylated Counts
#' @param training.chr.sel the chromsome that should be used for the training of the HMM
#' @param seqLengths a vector representing the sequence length for each chromosome of the assembly
#' @param knn a vector of values to be used with KNN
#' @param q use quantile based distance cutoff value, for each chromsome, 1=disabled
#' @param cutoff a vector of distances to be used for cutoff in KNN
#' @param num.cores to be used CPU cores. Multithreading handled by doSNOW cluster
#' @param markLowDensitySegm whenever regions with low density points 5\%<,>95\% should be marked seperately
#' @description Estimates PMD segments based on alphaValue estimation combinded with hidden markov model prediction. 
#' The method is adapted from MethylseekR(MSR) and extended with KNN datapoint selection. 
#' KNN will choose the neighrest neighorbing points based on the genomic distance. 
#' These points will be used for alpha value estimation, like in MSR.
#' Furthermore a distance cutoff is definded such that far away datapoints will not be considered for alpha estimation.
#' This method is applied for multiple values for K, and all estimated alpha values will be averaged.
#' The hidden markov model will than be trained and used with those alpha values.
#' @return GRange object with the PMD/nonPMD segments
#' @export
segmentPMDsKNN <-
function(m,
         training.chr.sel="chr2",
         pdfFilename = NULL,
         seqLengths = "hg19",
         knn=c(20,35,55),
         q=1,
         cutoff = c(50000,70000,90000),
         num.cores = 1,
         markLowDensitySegm = T
        ){ 

        message('Segmenting with these parameters\n\n',
                ' - KNN\t\t\t',paste(knn,collapse = ','),'\n',
                ' - cutoffs:\t\t',paste(cutoff,collapse = ','),'\n',
                ' - train.Chr:\t\t',training.chr.sel,'\n',
                ' - seqLengths:\t\t',paste(seqLengths,collapse = ","),'\n',
                ' - q:\t\t\t',q,'\n',
                ' - num.cores:\t\t',num.cores,'\n',
                ' - markLowDensitySegm:\t',markLowDensitySegm,'\n'
                )
       
  
         message('generating cluster of size',num.cores)
  
   
        cl <- makeCluster(num.cores)
        registerDoSNOW(cl)
        if(seqLengths == "hg19"){
                message('Using human genome assembly hg19')
                seqLengths = seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)
        }
        seqLengths = filterPresentChromosomesNSort(seqLengths,m)
        message('Calculating KNNs')
        knns.list = calculateKNNs(m,knn)
        message('Generating offsets from KNN')
        offsets = calculateKNNsOffsets(knns.list,cutoff,q)

        message('estimating alphaValues for each K')
        alphas = alphaValueEstimation(m,offsets)
        message('clearing up memory')
        remove(offsets)
        gc(verbose = F)
        message('combing alphaValues')
        
        df = do.call(cbind,alphas)
        comb_alpha = foreach(i = 1:nrow(df)) %dopar%{
                apply(sapply(df[i,],cbind),1,mean) 
        }
        #comb_alpha = lapply(comb_alpha,'/',length(knn))
        names(comb_alpha) = names(alphas[[1]])
        message('clearing up memory')
        remove(alphas)
        gc(verbose = F)
        
        hmm.model <- trainPMDHMM(comb_alpha, training.chr.sel, plot.distr=TRUE, pdfFilename)
        y.list <- PMDviterbiSegmentation(comb_alpha, hmm.model,T,knns.list[[1]])
        segments <- createGRangesObjectPMDSegmentation(m, y.list, seqLengths)
        message('shutting down cluster')
        stopCluster(cl)
        
        segments
         
 }


filterPresentChromosomesNSort = function(seqlengths,gr){
        seqlengths = seqlengths[names(seqlengths) %in% levels(seqnames(gr))]
        seqlengths[order(match(names(seqlengths),names(seqlengths(gr))))]
}
