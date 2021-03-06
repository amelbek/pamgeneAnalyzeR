% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{substractBackground}
\alias{substractBackground}
\title{Background substraction}
\usage{
substractBackground(rawData)
}
\arguments{
\item{rawData}{Raw data files as obtained from \code{\link{mergeExperiments}} function where peptides are in columns and samples in rows}
}
\value{
Return the normatized data frame where from each row the mean of all background points was substracted
}
\description{
Substract the backgroud for signal normalization
}
