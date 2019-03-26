#' Find blobs centers
#'
#' Find in the PamChip image the center of each peptide, left and right controls and 8 background points. Use the \code{\link{plot_centers}} function to check that the centers were detected correcly and adapt the threshold with the \code{thresh} parameter.
#'
#' @param image PamChip image.
#' @param thresh threshold the image to get only the pixels with the highest values. either "auto" or or percentage string (default "99.7%").
#'
#' @return Returns a data frame with the \code{x} and \code{y} coordinates of the centers, thier row and column number.
#' @export
find_centers <- function(image, thresh = "99.7%"){
  if (missing(image)){
    stop("use TIFF image as input")
  }
  if (class(image) != "matrix"){
    stop("use TIFF image as input")
  }
  #####################################################
  ## sub-fuction to get the centers of detected blobs
  #####################################################
  get.centers <- function(im,thr=thresh)
  {

    imh <- imager::imhessian(im) # get Hessian matrix
    dt <-  imh$xx*imh$yy - imh$xy^2 # get Derivative
    dt <- imager::threshold(dt, thr) # set threshold
    dt <- imager::label(dt)

    dt <- as.data.frame(dt)
    dt <- subset(dt, value>0)
    dt <- dplyr::group_by(dt, value)
    dt <- dplyr::summarise(dt, mx=mean(x),my=mean(y))
  }
  #####################################################


  # We use only left side control to detect the blobs
  subf1_left <- image[,0:(ncol(image)/4)]
  blobs <- imager::as.cimg(t(subf1_left))
  #plot(blobs)
  nblobs.denoised <- imager::isoblur(blobs,10, gaussian=TRUE)

  centers <- get.centers(nblobs.denoised)

  ## Use only the top left point as reference for other points
  ## The other points are at constant distance fro this point

  # Get peptides Centers
  xs <- seq(centers$mx[centers$my ==min(centers$my)]+86,centers$mx[centers$my ==min(centers$my)]+324, by=21.6)
  rows <- c(1:12)
  ys <- seq(min(centers$my)-65,min(centers$my)+175, by=21.6)
  columns <- c(1:12)

  centers.coords <- NULL

  for (x in 1:length(xs)){
    for (y in 1:length(ys)){
      centers.coords <- cbind(centers.coords, c(xs[x],ys[y], rows[x], columns[y]))
    }
  }

  colnames(centers.coords) <- paste("Peptide", 1:ncol(centers.coords), sep = "")

  # Get Controls centers
  #left
  centers.coords <- cbind(centers.coords, control_left1=c(centers$mx[centers$my ==min(centers$my)], min(centers$my), 1, -1))
  centers.coords <- cbind(centers.coords, control_left2=c(centers$mx[centers$my ==min(centers$my)]+23, min(centers$my)+44, 1, -1))
  centers.coords <- cbind(centers.coords, control_left3=c(centers$mx[centers$my ==min(centers$my)], min(centers$my)+44, 1, -1))
  centers.coords <- cbind(centers.coords, control_left4=c(centers$mx[centers$my ==min(centers$my)], min(centers$my)+87, 1, -1))
  #right
  centers.coords <- cbind(centers.coords, control_right1=c(centers$mx[centers$my ==min(centers$my)]+410, min(centers$my)+22, 1, -1))
  centers.coords <- cbind(centers.coords, control_right2=c(centers$mx[centers$my ==min(centers$my)]+388, min(centers$my)+108, 1, -1))
  centers.coords <- cbind(centers.coords, control_right3=c(centers$mx[centers$my ==min(centers$my)]+410, min(centers$my)+66, 1, -1))
  centers.coords <- cbind(centers.coords, control_right4=c(centers$mx[centers$my ==min(centers$my)]+410, min(centers$my)+108, 1, -1))

  ## backgrounds
  # left
  centers.coords <- cbind(centers.coords, left_background1=c(centers$mx[centers$my ==min(centers$my)]+30, min(centers$my),-3, -3))
  centers.coords <- cbind(centers.coords, left_background2=c(centers$mx[centers$my ==min(centers$my)]+30, min(centers$my)+88,-3, -3))

  # right
  centers.coords <- cbind(centers.coords, right_background1=c(centers$mx[centers$my ==min(centers$my)]+380, min(centers$my)+20,-4, -4))
  centers.coords <- cbind(centers.coords, right_background2=c(centers$mx[centers$my ==min(centers$my)]+380, min(centers$my)+66, -4, -4))

  #bottom
  centers.coords <- cbind(centers.coords, bottom_background1=c(centers$mx[centers$my ==min(centers$my)]+300, min(centers$my)+220, -2, -2))
  centers.coords <- cbind(centers.coords, bottom_background2=c(centers$mx[centers$my ==min(centers$my)]+100, min(centers$my)+220, -2, -2))

  #top
  centers.coords <- cbind(centers.coords, top_background1=c(centers$mx[centers$my ==min(centers$my)]+300, min(centers$my)-120, -1, -1))
  centers.coords <- cbind(centers.coords, top_background2=c(centers$mx[centers$my ==min(centers$my)]+100, min(centers$my)-120, -1, -1))

  row.names(centers.coords) <- c('y', 'x', 'Row', 'Col') # we inverse x and y since we use the transposed image first
  centers.coords <- as.data.frame(t(centers.coords))
  centers.coords$names <- row.names(centers.coords)
  row.names(centers.coords) <- NULL
  centers.coords <- centers.coords[c('names', 'x', 'y', 'Row', 'Col')]
  return(centers.coords)

}



###################################################################################

#' Plot centers
#'
#' Plot the PamChip image with the detecte centers from the \code{\link{find_centers}} function. This function is used to check that the centers are detected correcly.
#'
#' @param image PamChip image.
#' @param centers_coords data frame of coordinates of centers obtained from the \code{\link{find_centers}} function.
#'
#' @return Returns the PamChip image with the centers of each peptides, left ans right controls and 8 background points.

#' @export
plot_centers <- function(image, centers_coords){

  if (missing(image)){
    stop("use TIFF image as input")
  }
  if (missing(centers_coords)){
    stop("use find_centers output as input")
  }

  plot(imager::as.cimg(t(image)))

  peptides <-  centers_coords[grepl("Peptide", centers_coords$names),]
  points(peptides$y,peptides$x,col="red", pch = 3)

  controls <-  centers_coords[grepl("control", centers_coords$names),]
  points(controls$y,controls$x,col="green", pch = 3)

  backgrouds <- centers_coords[grepl("background", centers_coords$names),]
  points(backgrouds$y,backgrouds$x,col="blue", pch = 4)

  legend("bottomleft", legend=c("Peptides", "Controls", "Backgrounds"),
    col=c("red", "green", "blue"), pch=c(3,3,4), cex=0.8, box.lty=0)

  title(xlab="Y coordinates", ylab="X coodinates")
}


