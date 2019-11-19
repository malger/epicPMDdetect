writeSegmentation =function(granges,seg.file.name){
  red = c(255,95,95)
  green = c(92,229,85)
  orange = c(255, 165, 0)

  cchose = Vectorize(function(x) switch (x,
                                        "notPMD" = paste0(green,collapse=","),
                                        "PMD" = paste0(red,collapse=","),
                                        "Potential_PMD" = paste0(orange,collapse=",")
                                ))
  
  seg <- data.frame(seqnames=seqnames(granges),
                    starts=start(granges),
                    ends=end(granges),
                    compartment=(elementMetadata(granges)$type),
                    nCG = (elementMetadata(granges)$nCG),
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
