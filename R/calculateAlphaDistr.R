#' Estimates AlphaValues based on window offsets for a Chromosome'
#' @param m Granges of one particular chromosome. Metadata musst include Total,Methylated Counts
#' @param offsets list of relative offsets describing the datapoints to be used for alpha estimation
#' @return estimated alpha values
calculateAlphaDistr <-
function(M, T, offsets){
  
  binSize <- 0.1 # bin size to discretize alpha values
  alphas <- seq(binSize, 3, by=binSize) # discretized alphas
  
  logPs <- foreach(a = alphas) %do% rollapply((lbeta(M+a, T-M+a) - lbeta(a, a)), width=offsets,sum,align="left")
  
  logPs <- do.call("cbind", logPs)
  score <- apply(logPs, 1, function(x){m=max(x,na.rm = T); p=exp(x-m); p=p/sum(p); sum(p*alphas)})
  score
}
