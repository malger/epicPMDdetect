calculateAlphaDistr <-
function(M, T, num.cores, offsets){
  
  binSize <- 0.1 # bin size to discretize alpha values
  alphas <- seq(binSize, 3, by=binSize) # discretized alphas
  
  print('offsets finished')

  logPs <- mclapply(alphas, function(a){
    r = zoo::rollapply((lbeta(M+a, T-M+a) - lbeta(a, a)), width=offsets,sum,align="left")
    r
   # c(r,rep(r[length(r)],length(M)-length(r)))
  #  notnorm/sapply(offsets_bins,length)  # #logPS muss nach der anzahl zusammengefasster cpgs normalisiert werden 
  }, mc.cores=num.cores)

  logPs <- do.call("cbind", logPs)
  score <- apply(logPs, 1, function(x){m=max(x,na.rm = T); p=exp(x-m); p=p/sum(p); sum(p*alphas)})
  score
}
