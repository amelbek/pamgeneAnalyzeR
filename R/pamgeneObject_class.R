
#########################################################################################################
# Create a class pamgene to read fiels with signal and extraxt info about cycle and exposure time
#########################################################################################################
#' Pamgene object from the extracted signal .txt file
#'
#' Takes as input the .txt file and return an object containing the data for each experiment with the corresponding metadata: chipID, well, exposure time and cycle
#'
#' @param file the extracted signal .txt file name.
#' @return Takes as input the .txt file and return an object containing the data for each experiment with the corresponding metadata: chipID, well, exposure time and cycle
#'
#' @export
pamgeneObject <- function(file, data, chipID, well, exptime, cycle){

  if (grepl(".txt", file) == FALSE ){
    stop(" file must be text file")
  }

  attribs <- strsplit(basename(file), split = '_')
  d <- read.table(file, header = TRUE)

  d$ID <- as.character(d$ID)
  d$ID[is.na(d$ID)] <- as.character(d$names[is.na(d$ID)])

  z <- list(data= d, chipID = attribs[[1]][1],well= attribs[[1]][grep("^[W].*", attribs[[1]])], exposureTime= attribs[[1]][grep("^[T].*", attribs[[1]])], cycle=attribs[[1]][grep("^[P].*", attribs[[1]])] )
  class(z) <- append(class(z),"pamgeneObject")
  return(z)
}



#' merge experiments
#'
#' This function merge all experiments (from extracted signal files) in a directory in a data frame with peptides in columns and samples in rows. It also adds as columns the metadata chipID , well, cylce and exposure time.
#'
#' @param path the path for the directory where the .txt experiments file obtained from signal exptraction are.
#' @param pattern an optional regular expression. Only files names which match the regular expression will be used. if NULL all files in the directory will be used (default NULL).
#' @param value The value from the extracted signal files to be used. "median", "mean" or "sum" (default "median").
#' @return Return a data frame with with  chipID , well, cylce exposure time and petides in columns and samples in rows.
#'
#' @export
#'
mergeExperiments <- function(path, pattern = NULL, value = "median") {

  batch <- dir(path = path, pattern = pattern)
  length(batch)

  d <- NULL
  for (filename in batch){
    pam <- pamgeneObject(file = paste(path,filename, sep = '/'))

    tmp <- pam$data[,value]
    names(tmp) <- pam$data$ID
    d <- rbind(d, data.frame(chipID = pam$chipID, well = pam$well, cycle = pam$cycle, exposureTime = pam$exposureTime, t(tmp), check.names = FALSE))
  }

  d$cycle = factor(d$cycle,levels=as.character(unique(d$cycle))[gtools::mixedorder(as.character(unique(d$cycle)))])
  d$exposureTime = factor(d$exposureTime,levels=as.character(unique(d$exposureTime))[gtools::mixedorder(as.character(unique(d$exposureTime)))])

  return(d)
}




