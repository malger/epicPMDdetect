calculateKNNs = function(m,knn,cutoff=NULL,q=1){
  
  offsets = foreach(k = knn) %dopar% {
    message('calculating KNN for k=',k)
    chrs = unique(seqnames(m))
    k_offsets = foreach(chr.sel = chrs) %dopar% {
      chrindx=as.character(GenomicRanges::seqnames(m))==chr.sel;
      getKNNOffsets(GenomicRanges::start(m[chrindx]),k,q,cutoff)
    }
    names(k_offsets) = chrs
    k_offsets
  }
  names(offsets) = paste0('knn',knn)
  offsets
  
}
