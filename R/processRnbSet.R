segmentRnbSet = function(rnbobj,sample.col.name,outputFolder,num.cores=1,settings = NULL){
  data = readEPIC_RnbSet(rnbobj,sample.col.name)
  foreach (s = samples(data)) %do% {
    message('processing sample: ',s)
    m = buildMethGrangesFromRnbSet(data,s)
    segm = segmentPMDsKNN(m,num.cores)
    writeSegmentation(segm,file.path(outputFolder,'/',s,'epicPMDdetect'))
  }
}