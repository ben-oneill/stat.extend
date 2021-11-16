#' Logarithm of the subfactorial numbers
#'
#' \code{lsubfactorial} returns the logarithms of the subfactorial numbers.
#'
#' The subfactorial numbers count the number of derangements of a set of objects (permutations in which no element
#' appears in its original position).  This function computes the logarithms of the subfactorial numbers for a given
#' input vector specifying the numbers of interest.
#'
#' @param x A vector of non-negative integers
#' @return If the input is a vector of non-negative integers, the output will be a vector of the logarithms of
#' the corresponding subfactorial numbers.
#'
#' @examples 
#' # In the limit n! / !n goes to e
#' # so limit of differences of logs is 1
#' lfactorial(1000) - lsubfactorial(1000)
#' 
#' 

lsubfactorial <- function(x) {
  
  #Compute values of the subfactorials recursively
  MAX <- max(x, na.rm = TRUE)
  LL  <- rep(0, max(2,MAX+1))
  LL[2] <- -Inf
  if (MAX > 1) {
    for (k in 2:MAX) {
      LL[k+1] <- log(k-1) + matrixStats::logSumExp(c(LL[k], LL[k-1])) } }
  
  #Compute output
  OUT <- rep(0, length(x))
  for (i in seq_along(x)) {
    if (x[i] < 0  || is.na(x[i]) ) {
      OUT[i] <- NA } 
    else {
      OUT[i] <- LL[x[i]+1]}}
  
  #Return output
  OUT }