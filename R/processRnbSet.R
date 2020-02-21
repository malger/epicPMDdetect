#' Segments all samples of an Rnbeads Object
#' @param rnbobj the rnbeads object with the methylation data
#' @param outputFolder the folder to which the segmentation files should be written to
#' @param num.cores number of cores to use in parallel
#' @param settings pass settings to segmentation Function. Not yet implemented
#' @export
#' @import RnBeads
segmentRnbSet = function(rnbobj,outputFolder,num.cores=1,settings = NULL){
  for (i in 1:length(samples(data))) {
    s = samples(data)[i]
    message('processing sample: ',s)
    m = buildMethGrangesFromRnbSet(data,s)
    segm = segmentPMDsKNN(m,num.cores = num.cores)
    writeSegmentation(segm,file.path(outputFolder,paste0(s,'.segm')))
    writePMDs2bed(segm,file.path(outputFolder,paste0(s,'.pmds')))
    writeNpmds2bed(segm,file.path(outputFolder,paste0(s,'.oths')))
  }
}
