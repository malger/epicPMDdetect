% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processRnbSet.R
\name{segmentRnbSet}
\alias{segmentRnbSet}
\title{Segments all samples of an Rnbeads Object}
\usage{
segmentRnbSet(
  rnbobj,
  outputFolder,
  samples = NULL,
  num.cores = 1,
  settings = NULL
)
}
\arguments{
\item{rnbobj}{the rnbeads object with the methylation data}

\item{outputFolder}{the folder to which the segmentation files should be written to}

\item{samples}{optional param: if provided only those samples are processed}

\item{num.cores}{number of cores to use in parallel}

\item{settings}{pass settings to segmentation Function. Not yet implemented}
}
\description{
Segments all samples of an Rnbeads Object
}
