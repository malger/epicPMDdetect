% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readEPICidats.R
\name{readEPIC_RnbSet}
\alias{readEPIC_RnbSet}
\title{Reads an existing RnBeads Set object.}
\usage{
readEPIC_RnbSet(path, sample.col.name, preprocess = T)
}
\arguments{
\item{path}{to the rnbset object saved on disk}

\item{sample.col.name}{name of the column in the csv that contains the names of the samples}

\item{preprocess}{whenever rnbeads default preprocessing pipeline should be run. PP results are saved in pp_report/ relatived to workin dir}
}
\value{
rnbeadsSet object
}
\description{
Reads an existing RnBeads Set object.
}
