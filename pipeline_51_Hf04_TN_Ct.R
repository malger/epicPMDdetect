hmm.train.chr="chr2"
num.cores=6
library(msrEPIC)
filterPresentChromosomesNSort = function(seqlengths,gr){
  seqlengths = seqlengths[names(seqlengths) %in% levels(seqnames(gr))]
  
  seqlengths[order(match(names(seqlengths),names(seqlengths(gr))))]
}
rnb_data = readEPIC_RnbSet("~/epic_seg/data/51/EPIC/rnb.set/","SampleID",F) 
m = msrEPIC::buildMethGranges(rnb.obj = rnb_data,samples = "51_Hf04_TN_Ct")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
seqLengths = seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)
seqLengths = filterPresentChromosomesNSort(seqLengths,m)
hmm.model <- msrEPIC::trainPMDHMM(m, hmm.train.chr, num.cores, 
                         plot.distr=TRUE, pdfFilename=NULL,
                         method="knn",NULL,k=30,q=1,cutoff=NULL
                         )
HMMprediction <- PMDviterbiSegmentation(m,
                                        hmm.model, 
                                        num.cores,
                                        "knn",
                                        NULL,
                                        30,
                                        1,
                                        cutoff=NULL,
                                        markLowDensityPoints=F
                )
out = m
elementMetadata(out) = elementMetadata(m)[,c("T","M","ID","CGI Relation")]
out$hmmscore = rep(NA)
out$hmmpred = rep(NA)



for(chr in names(HMMprediction)){
  out[seqnames(out)==chr]$hmmscore = HMMprediction[[chr]]$score
  out[seqnames(out)==chr]$hmmpred = HMMprediction[[chr]]$s
}
write.csv(data.frame(out),file="msrKNN.30.51_Hf04_TN.hmmpred.csv")
