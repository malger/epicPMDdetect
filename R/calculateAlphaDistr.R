calculateAlphaDistr <-
function(M, T, offsets){
  
  binSize <- 0.1 # bin size to discretize alpha values
  alphas <- seq(binSize, 3, by=binSize) # discretized alphas
  
  print('offsets finished')

  logPs <- foreach(a = alphas) %dopar% zoo::rollapply((lbeta(M+a, T-M+a) - lbeta(a, a)), width=offsets,sum,align="left")
  
  logPs <- do.call("cbind", logPs)
  score <- apply(logPs, 1, function(x){m=max(x,na.rm = T); p=exp(x-m); p=p/sum(p); sum(p*alphas)})
  score
}
