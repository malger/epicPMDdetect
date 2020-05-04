default.options = list(
  training.chr.sel="chr2",
  pdfFilename = NULL,
  seqLengths = "hg19",
  knn=c(30,60,90),
  q=1,
  cutoff = c(50000,70000,90000),
  num.cores = 1,
  markLowDensitySegm = F
)

#' Segments all samples of an Rnbeads Object
#' @param rnbobj the rnbeads object with the methylation data
#' @param outputFolder the folder to which the segmentation files should be written to
#' @param samples optional param: if provided only those samples are processed
#' @param num.cores number of cores to use in parallel
#' @param settings pass settings to segmentation Function. Not yet implemented
#' @export
#' @import RnBeads,utils
segmentRnbSet = function(rnbobj,outputFolder,samples=NULL,num.cores=1,settings = NULL){
  if(is.null(samples))
    samples = samples(rnbobj)
  if (!all(samples %in% samples(data)))
    stop("some samples specified were not found in the rnbset. Aborting!")
  for (i in 1:length(samples(rnbobj))) {
    s = samples(rnbobj)[i]
    message('processing sample: ',s)
    m = buildMethGrangesFromRnbSet(rnbobj,s)
    if(is.null(settings)){
      segm = segmentPMDsKNN(m,num.cores = num.cores)
    } else {
      if(typeof(settings) != "list")
        stop("settings musst be a list with following named children: knn,cutoff,seqLengths")
      #TODO pass more settings
      args = c(list(m=m),modifyList(default.options,settings))
      args[["num.cores"]]=num.cores
      segm = do.call(segmentPMDsKNN,args)
    }
    writeSegmentation(segm,file.path(outputFolder,paste0(s,'.segm')))
    writePMDs2bed(segm,file.path(outputFolder,paste0(s,'.pmds')))
    writeNpmds2bed(segm,file.path(outputFolder,paste0(s,'.oths')))
  }
}
