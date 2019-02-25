
#' \eqn{Z' factor}
#'
#' The \eqn{Z' factor} is a measure of the assay quality, showing the separation between the distribution of the positive and negative controls. one use the background as negative control and the bright left and right points as positive controls.
#'
#' @param positiveControls vector of positive controls values. In The Pamgene context it is the left and right bright controls points.
#' @param negativeControls vector of negative controls values. In The Pamgene context it is background points.
#'
#' @return Return \eqn{Z' factor} value.
#' @export
#'
calculateZpfactor <- function(positiveControls, negativeControls){

  SDp = sd(as.numeric(as.vector(positiveControls)))
  SDn = sd(as.numeric(as.vector(negativeControls)))
  Mp = mean(as.numeric(as.vector(positiveControls)))
  Mn = mean(as.numeric(as.vector(negativeControls)))

  return(1-((3*(SDp+SDn))/abs(Mp-Mn)))
}

#################################################################################

#' \eqn{Z' factor} for all samples in a data frame
#'
#' Calculate the \eqn{Z' factor} for all PamChip samples (rows) in a data frame. The \eqn{Z' factor} is a measure of the assay quality, showing the separation between the distribution of the positive and negative controls. one use the background as negative control and the bright left and right points as positive controls.
#'
#'
#' @param rawData take raw data files as obtained from \code{\link{mergeExperiments}} function.
#' @return Returns the same input data frame with an added column of calculated \eqn{Z' factors} for each row.
#'

#' @export
ZpFactor_forAll <- function(rawData){

  n_controls <- c('right_background1', 'right_background2', 'left_background1', 'left_background2', 'top_background1', 'top_background2', 'bottom_background1', 'bottom_background2')
  p_controls <- c('control_right1', 'control_right2', 'control_right3', 'control_right4', 'control_left1', 'control_left2', 'control_left3', 'control_left4')

  rawData_Z <- data.frame(ZpFactor = apply(rawData, 1, function(x) {calculateZpfactor( positiveControls = x[p_controls], negativeControls = x[n_controls])} ), rawData )
  return(rawData_Z)

}

#################################################################################
#' Background substraction
#'
#' Substract the backgroud for signal normalization
#'
#' @param rawData Raw data files as obtained from \code{\link{mergeExperiments}} function where peptides are in columns and samples in rows
#' @return Return the normatized data frame where from each row the mean of all background points was substracted
#'

#' @export



substractBackground <- function(rawData) {

  pheno <- c("chipID", "well", "exposureTime", "cycle", "ZpFactor")
  n_controls <- c('right_background1', 'right_background2', 'left_background1', 'left_background2', 'top_background1', 'top_background2', 'bottom_background1', 'bottom_background2')
  sig <- rawData[,-which(names(rawData) %in% pheno)]
  if (ncol(sig) != 160) {
    stop("rawData input must contain all 140 petides + backgrounds and controls" )
  }
  n <- apply(sig,1, function(x) {(x- mean((x[n_controls])))} )

  out <- data.frame(rawData[,which(names(rawData) %in% pheno)], t(n))

  return(out)
}
