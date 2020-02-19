#' Save PMDs to bed file or RDS
#'
#' @param PMDs GRange segmentation
#' @param GRangesFilename path of the GRange output file
#' @param TableFilename path of the RDS output file
#' @description Saves PMD segments to disk. Use either GRangesFilename for bed file format or TableFilename for RDS format.
#' @export
savePMDSegments <-
function(PMDs, GRangesFilename = NULL, TableFilename = NULL){

  # save as GRanges object
  if(!is.null(GRangesFilename))
  saveRDS(PMDs, GRangesFilename)

  # save as tab-delimited table

  if(!is.null(TableFilename)){
    indx=values(PMDs)$type=="PMD"
    D=data.frame(chr=as.character(seqnames(PMDs))[indx], start=start(PMDs)[indx], end=end(PMDs)[indx])
    write.table(D, file=TableFilename, quote=FALSE, sep="\t", row.names=FALSE)
  }
  
}
