#' Write a PMD/nPMD Segmentation to a bed.gz file
#' @param granges representing the segmentation
#' @param file.name file path and name. Bed.gz will be concatenated
#' @param score.col.name an optional column that contains scores for each segment 
#' @export
writeSegmentation =function(granges,seg.file.name,score.col.name=NULL){
  red = c(239,141,137)
  blue = c(147,188,246)
  bluegrey =c(119,136,153)

  cchose = Vectorize(function(x) switch (as.character(x),
                                        "notPMD" = paste0(blue,collapse=","),
                                        "PMD" = paste0(red,collapse=","),
                                        "unkown" = paste(bluegrey,collapse = ",")
                                ))
  
  if(is.null(score.col.name)){
    score = rep(0)
  }else{
    score = elementMetadata(granges)[,score.col.name]
    
  }
    
  seg <- data.frame(seqnames=seqnames(granges),
                    starts=start(granges),
                    ends=end(granges),
                    compartment=(elementMetadata(granges)$type),
                    score = score,
                    strands=strand(granges),
                    rstarts=start(granges),
                    rends=end(granges),
                    colors = cchose(granges$type)
                    # cor = (elementMetadata(a_b_comp)$cor.matrix)
  )
  file.name = paste0(seg.file.name,".bed.gz")
  print("writing results to file")
  gzf = gzfile(file.name,open = "w")

  write.table(format(seg ,scientific=FALSE),
              file=gzf,
              quote=F,
              sep="\t", row.names=F, col.names=F
  )
  close(gzf)
}

#' Write a GRange object a bed.gz file
#' @param granges representing the segmentation
#' @param file.name file path and name. Bed.gz will be concatenated
writeGranges2bed = function(granges,seg.file.name){
  seg <- data.frame(seqnames=seqnames(granges),
                    starts=start(granges),
                    ends=end(granges)
  )
  file.name = paste0(seg.file.name,".bed.gz")
  print("writing results to file")
  gzf = gzfile(file.name,open = "w")
  
  write.table(format(seg ,scientific=FALSE),
              file=gzf,
              quote=F,
              sep="\t", row.names=F, col.names=F
  )
  close(gzf)
}

#' Write PMDs from a GRange Segmentation object a bed.gz file
#' @param granges representing the segmentation
#' @param file.name file path and name. Bed.gz will be concatenated
#' @export
writePMDs2bed = function(granges,seg.file.name){
  writeGranges2bed(granges[granges$type %in% c('likelyPMD','PMD')],seg.file.name)
}

#' Write nonPMDs from a GRange Segmentation object a bed.gz file
#' @param granges representing the segmentation
#' @param file.name file path and name. Bed.gz will be concatenated
#' @export
writeNpmds2bed = function(granges,seg.file.name){
  writeGranges2bed(granges[granges$type %in% c('likelynotPMD','notPMD')],seg.file.name)
}
