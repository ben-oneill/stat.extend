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


lsubfactorial <- function(x) {
  
  #Check input
  if (!is.numeric(x))                       stop('Error: Input x should be composed of non-negative integers')
  xx <- as.integer(x)
  if (!all(x == xx))                        stop('Error: Input x should be composed of non-negative integers')
  if (min(x) < 0)                           stop('Error: Input x should be composed of non-negative integers')
  
  #Compute values of the subfactorials recursively
  MAX <- max(xx)
  LL  <- rep(0, max(2,MAX+1))
  LL[2] <- -Inf
  if (MAX > 1) {
    for (k in 2:MAX) {
      LL[k+1] <- log(k-1) + matrixStats::logSumExp(c(LL[k], LL[k-1])) } }
  
  #Compute output
  OUT <- rep(0, length(xx))
  for (i in 1:length(xx)) {
    OUT[i] <- LL[xx[i]+1] }
  
  #Return output
  OUT }