#' @title  K-nearest neighbor density estimation using R
#' @description  the k-nearest neighbor density estimation for two-dimensional value using R, and this function needs R_package "FNN"
#' @param x data values which is two dimensional
#' @param k the order of the nearest neighbor
#' @param xrange,yrange  the range of each dimension value
#' @return estimated points' correspoding densities 
#' @import FNN
#' @examples
#' \dontrun{
#' x <- as.matrix(faithful)
#' xrange <- seq(from = 1, to = 6, by = 0.1) 
#' yrange <- seq(from = 40, to = 100, by = 0.5)
#' k <- 5 
#' fit <- knnde(x, k, xrange, yrange)
#' }
#' @export
knnde <- function(x, k, xrange, yrange){ 
  p <- ncol(x) 
  n <- nrow(x) 
  est_pt <- expand.grid(xrange, yrange) 
  distance <- FNN::knnx.dist(x, est_pt, k) 
  est <- matrix(k / (2 * n * distance[,k]), nrow = length(xrange)) 
  return(est) 
}
