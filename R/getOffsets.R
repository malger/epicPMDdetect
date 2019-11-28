getBinSizeOffsets = function(CpGbins,q){
  
  offsets = as.list(rep(0,length(M))) #as.list(-(1:(length(M))))
  j = 1 
  i = 1 
  while (j <= length(CpGbins)){
    # if(CpGbins[j]<10){
    #   offsets[i:(i+CpGbins[j])] = rep(list(0:10),CpGbins[j]+1)
    # }else{
    # offsets[i:(i+max(offsets_bins[[j]]))] = rep(offsets_bins[j],max(offsets_bins[[j]])+1)
    obsn = round(q*CpGbins[j]-1)
    offsets[i:(i+CpGbins[j])] = rep(list(0:obsn),CpGbins[j]+1)
    # }
    
    i=i+CpGbins[j]+1
    j=j+1
  }
  offsets
}

getKNNOffsets = function(positions,k,q=1,cutoff=NULL){
  nn = FNN::get.knnx(positions,positions,k,algorithm = 'kd_tree')
  offsets = apply(nn$nn.index,2,function(x) x-nn$nn.index[,1])
  out = split(offsets,1:nrow(offsets))
  if(!is.null(cutoff)){
    print(paste0("using cutoff of ",cutoff))
    dists = apply(nn$nn.dist,2,function(x) x-nn$nn.dist[,1])
    a = mclapply(1:nrow(dists),function(i)dists[i,]<cutoff)
    out = mclapply(1:length(out),function(i)out[[i]][a[[i]]])
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

