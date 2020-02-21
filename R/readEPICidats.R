#' Reads EPIC arrays idat files with RnBeads default pipeline
#' @param sample.sheet a csv describing the present samples.
#' @param idat.dir path where the idat files are located
#' @param sample.col.name name of the column in the csv that contains the names of the samples
#' @param preprocess whenever rnbeads default preprocessing pipeline should be run. PP results are saved in pp_report/ relatived to workin dir
#' @return rnbeadsSet object
#' @export
#' @import RnBeads
readEPIC_idats = function(sample.sheet,idat.dir,sample.col.name,preprocess=T){

  rnb.options(identifiers.column = sample.col.name)
  ds = c(idat.dir,sample.sheet )
  print("import data...")
  data = rnb.execute.import(ds,data.type="infinium.idat.dir")

  if(preprocess){
    rnb.run.preprocessing(rnb.obj,dir.reports=file.path(getwd(),'pp_report/'))$rnb.set
  } else {
    data
  }
}

#' Reads an existing RnBeads Set object.
#' @param path to the rnbset object saved on disk
#' @param sample.col.name name of the column in the csv that contains the names of the samples
#' @param preprocess whenever rnbeads default preprocessing pipeline should be run. PP results are saved in pp_report/ relatived to workin dir
#' @return rnbeadsSet object
#' @export
#' @import RnBeads
readEPIC_RnbSet = function(path,sample.col.name,preprocess=T){
  rnb.options(identifiers.column = sample.col.name)
  print("import data...")
  d = load.rnb.set(path)
 # print(paste("found following samples: ",samples(d)))
  d
}

#' Reads an Grange object saved in RDS format.
#' @param rdsfile to the rnbset object saved on disk
#' @return grange Object
#' @description The grange object musst contain (T)otal,(M)ethylated Read counts in the metadata. 
#' @export
readEPIC_RDSgranges = function(rdsfile){
  gr = readRDS(rdsfile)
  if(!c('M','T') %in% names(elementMetadata(gr))) {
    warning('Your Granges Object should contain an Value for (M)ethylation Bead Counts, and (T)otal Bead counts (U+M)')
    stop('Please provide T,M in your granges object. Aborting!')
  }
  gr = gr[!(is.na(gr$M))]
  #assure correct ordering
  gr <- gr[order(as.vector(seqnames(gr)), start(gr))]
  gr
  

}

#' Builds a Methylation GRange object from RnbSet that can be used for PMD detection.
#' @param rnb.obj a rnbeads object that represents an EPIC array experiment
#' @param samples the sample that will be used to create tbe Meth. Grange. If a vector is supplied the mean Meth over all specified samples is taken.
#' @return grange Object with Methylation and Total counts
#' @export
#' @import RnBeads
buildMethGrangesFromRnbSet = function(rnb.obj,samples=NULL){

  print("retrieving annotation...")
  gr = GRanges(annotation(rnb.obj,add.names = T))
  #filter existing chromosomes

  print("getting T,M Matrices...")
  unmethylated = U(rnb.obj,row.names = T)
  methylated = M(rnb.obj,row.names = T)
  betavals = meth(rnb.obj,row.names = T)

  C = unmethylated+methylated
  M = matrixcalc::hadamard.prod(betavals,C)
  remove(unmethylated,methylated,betavals)

  if(is.null(samples)){
    print("WARNING: Use Mean of ALL Samples!")

    M = apply(M,1,mean)
    C = apply(C,1,mean)
    values(gr) = cbind(data.frame(T=C,M),values(gr))
  }
  else if(length(samples) == 1){
    values(gr) = cbind(data.frame(T=C[,samples],M = M[,samples]),values(gr))
  } else {
    M = apply(M[,samples],1,mean)
    C = apply(C[,samples],1,mean)
    values(gr) = cbind(data.frame(T=C,M),values(gr))
    print("WARNING: Use Mean over specifed Samples!")
  }
  remove(M,C)
  #remove rows with NAs
  gr = gr[!(is.na(gr$M))]
  #assure correct ordering
  gr <- gr[order(as.vector(seqnames(gr)), start(gr))]
  gr
}

# calcDistanceNCpG = function(gr){
#   lapply(unique(seqnames(granges)),FUN=function(s) {
#     gr.chr = granges[seqnames(granges)==s]
#     c((start(gr.chr[-1])-end(gr.chr[-length(gr.chr)])),NA)
#   })
# }

# ADD R.utils to dependcies for this!
# filterOpenSeaProbes = function(gr){
#   printf("#CpGs before filtering OpenSea probes : %d\n",length(gr))
#   nonopenSeaIds = gr$`CGI Relation`  %in%  "Open Sea"
#   gr = gr[nonopenSeaIds,]
#   printf("#CpGs after filtering OpenSea probes: %d\n",length(gr))
#   gr
# }
