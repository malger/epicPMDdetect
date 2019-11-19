

getDynamicCpGBins =function(cpgs.gr,tilewidth = 250000){
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(txdb),
                                    tilewidth=tilewidth,
                                    cut.last.tile.in.chrom=TRUE
  )
  chromosomes = intersect(unique(seqnames(cpgs.gr)),unique(seqnames(bins)))
  o = lapply(chromosomes,FUN=function(x){
    print(x)
    s.bins = bins[seqnames(bins)==x]
    s.cpgs = cpgs.gr[seqnames(cpgs.gr)==x]
    print("ok")
    cpgs.ov = findOverlaps(s.bins,s.cpgs)
    nonEmptybins = s.bins[queryHits(cpgs.ov)]
    
    nonEmptybins$cpgs = names(s.cpgs[subjectHits(cpgs.ov)])
    nonEmptybins = nonEmptybins[!duplicated(nonEmptybins$cpgs)]
    uniqueNEbins = unique(nonEmptybins)
    countOverlaps(uniqueNEbins,nonEmptybins,type="equal")
  })
  names(o) = chromosomes
  o  
  
}