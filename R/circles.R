


# This function plot the circles from wich the signal will be taken
# can be used to test different radius
# take as input the centers found from find_centers

#' Plot circles on PamChip image
#'
#' Plot the circles on the PamChip image of the areas from wich the signal will be taken with the \code{\link{extract_signal}} function. this function is used to optimize the radius to be used and check that the circles are inside the blobs.
#'
#' @param img PamChip image.
#' @param centers.coords data frame of coordinates of centers obtained from the \code{link{find_centers}} function.
#' @param radius the number of pixels to take from the center to build the circle (default 5).
#' @param parallel Logical. If TRUE parallelize to process to fasten the output (default TRUE).
#'
#' @return Returns the PamChip image with the circles of areas of each peptides, left ans right controls and 8 background points.


#' @export
#'
plot_circles <- function(centers.coords, img, radius = 5, parallel=TRUE){

  ############################################
  # Function for one circle
  ############################################
  one_circle <- function(img, x, y, radius) {
    point <- NULL
    for( i in 1:ncol(img)){
      for(j in 1:nrow(img)){
        if((x-j)^2+(y-i)^2 < radius^2){
          point <- rbind(point, c(j,i))
        }
      }
    }
    return(point)
  }
  ############################################

  if (parallel==FALSE){
    all_points <- NULL
    for(i in 1:nrow(centers.coords)){
      point <- one_circle(img = img, x=centers.coords[i,'x'], y=centers.coords[i,'y'], radius = radius )
      all_points <- rbind(all_points, point)
    }
    colnames(all_points) <- c('j', 'i')
    all_points <- as.data.frame(all_points)
    plot(imager::as.cimg(t(img)))
    points(all_points$i,all_points$j,col="red", pch = '.')

  } else if (parallel==TRUE){
    cores=parallel::detectCores()
    cl <- parallel::makeCluster(cores-1) #not to overload your computer
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`
    all_points <- foreach::foreach(i=1:nrow(centers.coords), .combine=rbind) %dopar% {
      point <- one_circle(img = img, x=centers.coords[i,'x'], y=centers.coords[i,'y'], radius = radius )
      point
    }
    parallel::stopCluster(cl)
    colnames(all_points) <- c('j', 'i')
    all_points <- as.data.frame(all_points)
    plot(imager::as.cimg(t(img)))
    points(all_points$i,all_points$j,col="red", pch = '.')

  }

}




##########################################
# function to extract signal from whole image
##########################################

#' Extract raw signal from PamChip image
#'
#' function to extract raw signal of each peptide, left and right controls and 8 background points.
#'
#' @param img Pamgene image.
#' @param centers.coords data frame of coordinates of centers obtained from the \code{find_centers} function.
#' @param radius the number of pixels to take from the center to build the circle (default 5).
#' @param parallel Parallelize to process to fasten the output.
#'
#' @return Returns a data frame with the basic statistics values of all pixels in the circle (mean, median, sum ad sd) of each peptides, left ans right controls and 8 background points.

#' @export
#'
extract_signal <- function(centers.coords, img, radius = 5, parallel=FALSE){
  ###########################################
  ## States for one circle
  ###########################################
  circle.stats<-function(img,x,y,radius) {


    #coor <- reshape2::melt(img) #retrieve coodinates from the image
    #print(head(coor))

    ### find the points that are inside the circle
    #insiders <- apply(coor, 1, function(co) (x-co[1])^2+(y-co[2])^2 < radius^2)
    #coords.insiders <- coor[insiders,]

    coords.insiders <- NULL
    for( i in 1:ncol(img)){
      for(j in 1:nrow(img)){
        if((x-j)^2+(y-i)^2 < radius^2){
          coords.insiders <- rbind(coords.insiders, img[j,i])
        }
      }
    }

    coords.insiders <- as.data.frame(coords.insiders)
    colnames(coords.insiders) <- c('value')
    #print(coords.insiders)
    #print(dim(coords.insiders))
    median <- median(coords.insiders$value)
    mean <- mean(coords.insiders$value)
    max <- max(coords.insiders$value)
    min <- min(coords.insiders$value)
    sum <- sum(coords.insiders$value)
    var <- var(coords.insiders$value)
    #return(img)
    return(as.data.frame(cbind(mean,median,max,min,sum,var)))

  }
  ############################################################
  if(parallel==FALSE) {
    signal <- NULL

    for(row in 1:nrow(centers.coords)){
      temp <- cbind((centers.coords[row,]), circle.stats(x=centers.coords[row,'x'], y=centers.coords[row,'y'], radius= radius, img = img))
      signal <- rbind(signal, temp)

    }

    #parralelize for one image
    } else if (parallel==TRUE) {
      cores=parallel::detectCores()
      cl <- parallel::makeCluster(cores-1) #not to overload the computer
      doParallel::registerDoParallel(cl)
      `%dopar%` <- foreach::`%dopar%`
      signal <- foreach::foreach(row=1:nrow(centers.coords), .combine=rbind) %dopar% {
      temp = cbind((centers.coords[row,]), circle.stats(x=centers.coords[row,'x'], y=centers.coords[row,'y'], radius= radius, img = img))
      temp
      }
      parallel::stopCluster(cl)

    }
  return(signal)
}




#####################################################
## Function to plot the signal as the original image
## to check that it looks like the original image
#####################################################

# Take as input the extract_signal function output

#' Plot extracted signal as a PamChip image
#'
#' Plot the extracted signal from \code{\link{extract_signal}} function in the form of the original PamChip image. Used check that the signal was extracted correctly. The resulting plot should looks like the original PamChip image.
#'
#' @param signal data frame from the \code{\link{extract_signal}} function output.
#'
#' @return Returns the Pamgene image with the circles of areas of each peptides, left ans right controls and 8 background points. The resulting plot should looks like the original PanChip image.


#' @export
#'
plot_image_from_signal <- function(signal){

  g= min(signal$median)/max(signal$median)
  backround_color <- grey(g)

  ggplot2::ggplot(signal, ggplot2::aes(x=y, y=x, color=median) )+ ggplot2::scale_y_reverse() +ggplot2::geom_point(ggplot2::aes(size = 15))+
    ggplot2::scale_colour_gradientn(colours=c(backround_color,'white'))+
    ggplot2::theme(panel.background=ggplot2::element_rect(fill=backround_color), panel.grid=ggplot2::element_blank())+
    ggplot2::theme(axis.line=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
      axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), legend.position="none",
      panel.border=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
      plot.background=ggplot2::element_blank())

}



