getKNNOffsets = function(positions,k,q=1,cutoff=NULL){
  nn = FNN::get.knnx(positions,positions,k,algorithm = 'kd_tree')
  offsets = apply(nn$nn.index,2,function(x) x-nn$nn.index[,1])
  out = split(offsets,1:nrow(offsets))
  if(!is.null(cutoff)){
    print(paste0("using cutoff of ",cutoff))
    dists = apply(nn$nn.dist,2,function(x) x-nn$nn.dist[,1])
    a = mclapply(1:nrow(dists),function(i)dists[i,]<cutoff)
    out = mclapply(1:length(out),function(i)out[[i]][a[[i]]])
    
    
    #fix values that have under 5 datapoints
    useAdjNBs = function(x)unique(unlist(lapply(-2:2,function(i){
      x = ifelse(x<3,x+3,x)
      x = ifelse(x>length(out)-3,x-3,x)
      sapply(out[[(x+i)]],"+",i)
    })))
    ids = which(sapply(out,length)<5)
    for(i in ids){
      out[[i]] = useAdjNBs(i)
    }
    
    return(out)
  }
  
  if(q<1){
    print(paste0("using quantile of ",q))
    dists = apply(nn$nn.dist,2,function(x) x-nn$nn.dist[,1])
    distsl = split(dists[,-1],1:nrow(dists))
    maxdist = mclapply(distsl,quantile,probs=q,mc.cores=4)
    a = mclapply(1:nrow(dists),function(i)dists[i,]<maxdist[[i]])
    out = mclapply(1:length(out),function(i)out[[i]][a[[i]]])
  } else {
    out
  }
  
} 
getHighVariableDistNNs = function(nn) {
  nn.dist_var = apply(nn$nn.dist[,-1],1,mean)
  getOutliers <- function(x, removeNA = TRUE) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = removeNA)
    caps <- quantile(x, probs=c(.05, .95), na.rm = removeNA)
    H <- 1.5 * IQR(x, na.rm = removeNA)
    (x > (qnt[2] + H)) | (x<qnt[1] - H)
  }
  getOutliers(nn.dist_var)
}
