% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pamgeneObject_class.R
\name{mergeExperiments}
\alias{mergeExperiments}
\title{merge experiments}
\usage{
mergeExperiments(path, pattern = NULL, value = "median")
}
\arguments{
\item{path}{the path for the directory where the .txt experiments file obtained from signal exptraction are.}

\item{pattern}{an optional regular expression. Only files names which match the regular expression will be used. if NULL all files in the directory will be used (default NULL).}

\item{value}{The value from the extracted signal files to be used. "median", "mean" or "sum" (default "median").}
}
\value{
Return a data frame with with  chipID , well, cylce exposure time and petides in columns and samples in rows.
}
\description{
This function merge all experiments (from extracted signal files) in a directory in a data frame with peptides in columns and samples in rows. It also adds as columns the metadata chipID , well, cylce and exposure time.
}
