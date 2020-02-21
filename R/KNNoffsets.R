#' Calculates the offsets for each datapoint used to estimate it's alpha value
#' @param knns.list list(knnValue->list(chr)->KNNs)
#' @param cutoff distance cutoff at which datapoints should be discarded
#' @param q quantile based distance cutoff at which datapoints should be discarded, q=1 -> disabled
#' @return list(knnValue->list(chrs->vector(list(dp->offsets))))
#' @export
calculateKNNsOffsets = function(knns.list,cutoff=NULL,q=1){
  
  if(is.null(cutoff)) {
    cutoff = rep(list(NULL),length(knns.list))
  }
  
  chrs = names(knns.list[[1]])
  offsets = foreach(knns = knns.list,co = cutoff,.export='getKNNOffsets',.noexport="m") %:% 
    foreach(knns.chr = knns,.final = function(x) setNames(x, chrs)) %dopar% {
      getKNNOffsets(knns.chr,co,q)
    }
  
  names(offsets) = names(knns.list)
  offsets
  
}
#' Calculates the KNNs for each datapoint based on the genomic distance
#' @param m GenomicRange containg EPIC datapoints with Methylation, Total Value Counts
#' @param ks vector of values for K
#' @return list(knnValue->list(chrs->FNN knn-objects))
#' consider reading FNN::get.knnx documentation
calculateKNNs = function(m,ks){
  
  chrvec = seqnames(m)
  chrs = unique(chrvec)
  positions = start(m)
  
  nns = foreach(k = ks,.packages=c('GenomicRanges','FNN','mmap'),.export = 'as.boolean',.noexport = 'm') %:% 
    foreach(chr.sel = chrs,.final = function(x) setNames(x, chrs)) %dopar% {
      chrindx = as.boolean(chrvec==chr.sel)
      get.knnx(positions[chrindx],positions[chrindx],k,algorithm = 'kd_tree')
    }
  names(nns) = paste0('knn',ks)
  nns
  
}

getKNNOffsets = function(nn,cutoff =NULL,q=1){
  offsets = apply(nn$nn.index,2,function(x) x-nn$nn.index[,1])
  out = split(offsets,1:nrow(offsets))
  if(!is.null(cutoff)){
    print(paste0("using cutoff of ",cutoff))
    dists = apply(nn$nn.dist,2,function(x) x-nn$nn.dist[,1])
    a = lapply(1:nrow(dists),function(i)dists[i,]<cutoff)
    out = lapply(1:length(out),function(i)out[[i]][a[[i]]])
    
    
    # #fix values that have under 5 datapoints
    # useAdjNBs = function(x)unique(unlist(lapply(-2:2,function(i){
    #   x = ifelse(x<3,x+3,x)
    #   x = ifelse(x>length(out)-3,x-3,x)
    #   sapply(out[[(x+i)]],"+",i)
    # })))
    # ids = which(sapply(out,length)<5)
    # for(i in ids){
    #   out[[i]] = useAdjNBs(i)
    # }
    
    return(out)
  }
  
  if(q<1){
    print(paste0("using quantile of ",q))
    dists = apply(nn$nn.dist,2,function(x) x-nn$nn.dist[,1])
    distsl = split(dists[,-1],1:nrow(dists))
    maxdist = lapply(distsl,quantile,probs=q)
    a = lapply(1:nrow(dists),function(i)dists[i,]<maxdist[[i]])
    out = lapply(1:length(out),function(i)out[[i]][a[[i]]])
  } else {
    out
  }
  
}
#' Filters for data points that have a high variance in distance to their NNs
#' @param nn neighrest neighbor FNN::get.knnx object
#' @return indicies of the detected data points
#' @description Filters for datapoints that are out of the .25-.75 quantile in regards to 
#' their variance in distance
#' @export
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
