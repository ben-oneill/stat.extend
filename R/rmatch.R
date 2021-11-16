#' Random generation from the generalised matching distribution
#'
#' \code{rmatch} returns the a vector of pseudo-random values from the distribution.
#'
#' This function generates pseudo-random values from the generalised matching distribution.  Further details on 
#' the distribution can be found in the following paper:
#'
#' O'Neill, B. (2021) A generalised matching distribution for the problem of coincidences.
#'
#' @param n A non-negative integer specifying the number of rdandom values to produce
#' @param size The size parameter for the generalised matching distribution (number of objects to match)
#' @param prob The probability parameter for the generalised matching distribution (probability of known match)
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will
#' be a vector of pseudo-random numbers from the distribution

rmatch <- function(n, size, prob = 0) {
  
  #Check inputs
  if (!is.numeric(n))                        stop('Error: Input n should be numeric')
  if (!is.numeric(size))                     stop('Error: Size parameter should be a non-negative integer')
  if (length(n) != 1)                        stop('Error: Input n should be a single non-negative integer')
  if (length(size) != 1)                     stop('Error: Size parameter should be a single non-negative integer')
  if (size < Inf) { nn <- as.integer(size) } else { nn <- Inf }
  if (nn != size)                            stop('Error: Size parameter should be a non-negative integer')
  if (n < 0)                                 stop('Error: Input n should be a single non-negative integer')
  
  #Deal with trivial case where n = 0
  if (n == 0) { return(numeric(0)) }
  
  #Generate output vector
  OUT   <- rep(0, n)
  KNOWN <- rbinom(n, size = size, prob = prob)
  for (i in 1:n) {
    nnn    <- nn-KNOWN[i]
    PERM   <- sample.int(nnn, size = nnn, replace = FALSE)
    OUT[i] <- KNOWN[i] + sum(PERM == 1:nnn) }
  
  #Return output
  OUT }
