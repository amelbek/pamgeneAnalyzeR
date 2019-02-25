
#################################
# Function to register one image
# and extract signal from it
#################################



#' Align one image on a target image an extract the signal form it
#'
#' Uses the \code{\link{[RNiftyReg]}} package to align (register) the source image on the chosen target then extract the signal from it. the goal is to use one target image for which one extracted the coordinates of centers with the \code{\link{find_centers}} function and use this image as a reference on which one will align the pther images and extract the signal rom the same coordinates
#'
#' @param source the PamChip image to be aligned.
#' @param target The PamChip image to align on.
#' @param nLevels \code{\link[RNiftyReg]{niftyreg.nonlinear}} parameter (default 4).
#' @param maxIterations \code{\link[RNiftyReg]{niftyreg.nonlinear}} parameter (default 5).
#' @param useBlockPercentage \code{\link[RNiftyReg]{niftyreg.nonlinear}} (parameter default 50).
#' @param extract.signal if TRUE the basic statistics values of each peptides, left ans right controls and 8 background points will be extracted. if FALSE only the aligned (registered) image will be returned.
#' @param centers.coords data frame of coordinates of centers of the target image obtained from the \code{\link{find_centers}} function.
#' @param radius the number of pixels to take from the center to build the circle (default 5).
#' @param parallel.sig.extract if TRUE the process will be parallelized (default TRUE).
#' @param pamgene.layout.file STK or PTK PamChip layout.
#'
#' @return Returns a list with the registered image an the a data frame of extracted signal from it if the \code{extract.signal} parameter is TRUE. Only the registered image otherwise.
#'
#' @export
regiter_one_image <- function(source, target, nLevels=4, maxIterations=5,useBlockPercentage=50, extract.signal=TRUE, centers.coords, radius = 5, parallel.sig.extract = TRUE, pamgene.layout.file){

  if (missing(target))
    stop("specify a target image for registration")

  if (extract.signal==TRUE){

    transform <- RNiftyReg::niftyreg(source, target, nLevels=nLevels, maxIterations=maxIterations,useBlockPercentage=useBlockPercentage)
    resultImage <- RNiftyReg::applyTransform(RNiftyReg::forward(transform), RNiftyReg::retrieveNifti(source))

    signal <- extract_signal(centers.coords = centers.coords, img = resultImage, radius = radius, parallel = parallel.sig.extract)

    if (missing(pamgene.layout.file)){
      return(list(resultImage = resultImage, signal = signal))
    } else if (!missing(pamgene.layout.file)){

      signal2 <- merge(pamgene.layout.file, signal, by = c('Row', 'Col'), all = TRUE)
      return(list(resultImage = resultImage, signal = signal2))

    }

  } else if (extract.signal==FALSE){

    transform <- RNiftyReg::niftyreg(source, target, nLevels=nLevels, maxIterations=maxIterations,useBlockPercentage=useBlockPercentage)
    resultImage <- RNiftyReg::applyTransform(RNiftyReg::forward(transform), RNiftyReg::retrieveNifti(source))
    return(resultImage)
  }
}


#################################
# Function plot and check the registration
#################################

#' Plot PamChip images alignment
#'
#' Plot the registration the source on the target image to check that the alignment was done correctly.
#'
#' @param source the Pamgene image to be aligned.
#' @param target The Pamgene image to align on.
#' @param nLevels \code{\link[RNiftyReg]{niftyreg.nonlinear}} parameter (default 4).
#' @param maxIterations maxIterations \code{\link[RNiftyReg]{niftyreg.nonlinear}} parameter (default 5).
#' @param useBlockPercentage \code{\link[RNiftyReg]{niftyreg.nonlinear}} parameter (default 50).
#'
#' @return Returns the target image with on top the edges of the source image. The red edges should not be shifted.



#' @export
#'
check_registration <- function(source, target, nLevels=4, maxIterations=5,useBlockPercentage=50){
  transform <- RNiftyReg::niftyreg(source, target, scope='rigid', nLevels=nLevels, maxIterations=maxIterations,useBlockPercentage=useBlockPercentage)

  kernel <- mmand::shapeKernel(c(3,3))
  gradient <- mmand::dilate(transform$image,kernel) - mmand::erode(transform$image, kernel)
  mmand::display(target)
  mmand::display(mmand::threshold(gradient,method = 'kmeans'), add=T, col='red')

}




#################################
# Function to register all images in a given image
# and extract signal from them
#################################
#' @title Align all images in a directory on a target image and extract summary  statistics from each aligned image.

#' @usage Register (align) all the images in a directory on the chosen target.
#'
#' @param inDirectory The directory of the Pamgene images to be aligned.
#' @param outDirectory The output directoy to put the extracted signal and registered image on.
#' @param target The Pamgene image to align on.
#' @param nLevels RNiftyReg parameter.
#' @param maxIterations RNiftyReg parameter.
#' @param useBlockPercentage RNiftyReg parameter.
#' @param extract.signal if TRUE the basic statistics values of each peptides, left ans right controls and 8 background points will be extracted. if FALSE only the aligned (registered) image will be returned.
#' @param centers.coords data frame of coordinates of the target image of centers obtained from the \code{find_centers} function.
#' @param radius the number of pixels to take from the center to build the circle (default 5).
#' @param parallel.sig.extract if TRUE the signal extraction process will be parallelized (default FALSE).
#' @param pamgene.layout.file merge the Pamgene STK or PTK layout with the extracred values.
#' @param parallel.registration if TRUE the registration process will be parallelized.
#'
#' @return Returns for each image a .txt file with the extracted summary statistics for each peptide, left and right controls and 8 background points.
#'
#' @export
registerAllImages <- function(inDirectory, outDirectory, target, nLevels=4, maxIterations=5,useBlockPercentage=50, centers.coords, radius = 5, parallel.sig.extract = FALSE, pamgene.layout.file, parallel.registration = TRUE){
  dir.create(outDirectory)

  #################################################
  regiter_one_image <- function(source, target, nLevels, maxIterations,useBlockPercentage, extract.signal, centers.coords, radius, parallel.sig.extract, pamgene.layout.file){

    if (missing(target))
      stop("specify a target image for registration")

    if (extract.signal==TRUE){

      transform <- RNiftyReg::niftyreg(source, target, nLevels=nLevels, maxIterations=maxIterations,useBlockPercentage=useBlockPercentage)
      resultImage <- RNiftyReg::applyTransform(RNiftyReg::forward(transform), RNiftyReg::retrieveNifti(source))

      signal <- extract_signal(centers.coords = centers.coords, img = resultImage, radius = radius, parallel = parallel.sig.extract)

      if (missing(pamgene.layout.file)){
        return(list(resultImage = resultImage, signal = signal))
      } else if (!missing(pamgene.layout.file)){

        signal2 <- merge(pamgene.layout.file, signal, by = c('Row', 'Col'), all = TRUE)
        return(list(resultImage = resultImage, signal = signal2))

      }

    } else if (extract.signal==FALSE){

      transform <- RNiftyReg::niftyreg(source, target, nLevels=nLevels, maxIterations=maxIterations,useBlockPercentage=useBlockPercentage)
      resultImage <- RNiftyReg::applyTransform(RNiftyReg::forward(transform), RNiftyReg::retrieveNifti(source))
      return(resultImage)
    }
  }
  ####################################################

  if (missing(inDirectory)) {
    stop("specify images Directory path/name")
  }
  if (missing(inDirectory)) {
    stop("specify output Directory path/name")
  }
  if (missing(target)) {
    stop("specify a target image for registration")
  }

  ifelse(!dir.exists(file.path(outDirectory)), dir.create(file.path(outDirectory)), FALSE)
  ifelse(!dir.exists(file.path(outDirectory)), dir.create(file.path(outDirectory)), FALSE)

  file.names <- dir(inDirectory, pattern =".tif", full.names=TRUE) # add or .tiff

  if (parallel.registration == FALSE){

    for(i in 1:(length(file.names))){
      src <- tiff::readTIFF(source=file.names[i])
      reg <- regiter_one_image(source = src, target = target, nLevels=nLevels, maxIterations=maxIterations,useBlockPercentage=useBlockPercentage, extract.signal=TRUE, centers.coords=centers.coords, radius = radius, parallel.sig.extract = parallel.sig.extract, pamgene.layout.file = pamgene.layout.file)

      base_file_name <- sub("^([^.]*).*", "\\1", basename(file.names[i]))
      write.table(reg$signal, file = paste(outDirectory,'/' ,base_file_name, '.txt', sep = ''), sep = '\t',row.names = FALSE )


    }

  } else if (parallel.registration == TRUE){

    cores=parallel::detectCores()
    cl <- parallel::makeCluster(cores-1) #not to overload your computer
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`

    signal <- foreach::foreach(i=1:length(file.names), .export = "extract_signal") %dopar% {

      src <- tiff::readTIFF(source=file.names[i])
      reg <- regiter_one_image(source = src, target = target, nLevels=nLevels, maxIterations=maxIterations,useBlockPercentage=useBlockPercentage, extract.signal=TRUE, centers.coords=centers.coords, radius = radius, parallel.sig.extract = FALSE, pamgene.layout.file = pamgene.layout.file)

      base_file_name <- sub("^([^.]*).*", "\\1", basename(file.names[i]))
      write.table(reg$signal, file = paste(outDirectory,'/' ,base_file_name, '.txt', sep = ''), sep = '\t',row.names = FALSE )

    }
    parallel::stopCluster(cl)

  }
}


