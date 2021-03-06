% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trainPMDHMM.R
\name{trainPMDHMM}
\alias{trainPMDHMM}
\title{Trains and fits the Hidden Markov Model used for PMD assignment.}
\usage{
trainPMDHMM(alphas, chr.sel, plot.distr = TRUE, pdfFilename)
}
\arguments{
\item{alphas}{estimated alpha values. List of chromosomes->alphaValues}

\item{chr.sel}{the chromosome used for HMM training}

\item{plot.distr}{should the the postirior distributions for the HMM-states be plotted}

\item{pdfFilename}{save the former plot to pdf, default FALSE}
}
\value{
fitted Hidden Markov Model
}
\description{
Trains and fits the Hidden Markov Model used for PMD assignment. Method adapted from MethylseekR.
}
