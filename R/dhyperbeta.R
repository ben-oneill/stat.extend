#' Density function of the Gaussian hypergeometric beta distribution
#'
#' \code{dhyperbeta} returns the density or log-density values for the arguments.
#'
#' This function computes densities or log-densities from the density function of the directional form of the Gaussian 
#' hypergeometric beta distribution.  Further details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2025) Directional Gaussian hypergeometric beta distributions and their uses in contaminated binary sampling.
#'
#' @usage \code{dhyperbeta(x, shape1, shape2, intensity, proportion, right = FALSE, log = FALSE)}
#' @param x A vector of numeric values to be used as arguments for the density function
#' @param shape1 The shape1 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param shape2 The shape2 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param intensity The push-intensity parameter for the Gaussian hypergeometric beta distribution (non-negative numeric value)
#' @param proportion The push-proportion parameter for the Gaussian hypergeometric beta distribution (numeric value between zero and one)
#' @param right Logical direction value; \code{TRUE} uses the right-pushed Gaussian hypergeometric beta distribution; \code{FALSE} uses the left-pushed Gaussian hypergeometric beta distribution
#' @param log A logical value specifying whether results should be returned as log-densities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range and arguments are numeric)
#' then the output will be a vector of densities/log-densities corresponding to the vector argument x

dhyperbeta <- function(x, shape1, shape2, intensity, proportion, right = FALSE, log = FALSE) {

  #Check that inputs are appropriate type
  if (!is.numeric(x))                       stop('Error: Argument x is not numeric')
  if (!is.numeric(shape1))                  stop('Error: Shape1 parameter is not numeric')
  if (!is.numeric(shape2))                  stop('Error: Shape2 parameter is not numeric')
  if (!is.numeric(intensity))               stop('Error: Push-intensity parameter is not numeric')
  if (!is.numeric(proportion))              stop('Error: Push-proportion parameter is not numeric')
  if (!is.logical(right))                   stop('Error: Direction input (right) is not a logical value')
  if (!is.logical(log))                     stop('Error: log option is not a logical value')
  
  #Check that parameters and options are atomic
  if (length(shape1) != 1)                  stop('Error: Shape1 parameter should be a single number')
  if (length(shape2) != 1)                  stop('Error: Shape2 parameter should be a single number')
  if (length(intensity) != 1)               stop('Error: Push-intensity parameter should be a single number')
  if (length(proportion) != 1)              stop('Error: Push-proportion parameter should be a single number')
  if (length(right) != 1)                   stop('Error: Direction input (right) should be a single logical value')
  if (length(log) != 1)                     stop('Error: log option should be a single logical value')

  #Check that parameters are in allowable range
  if (min(shape1) <= 0)                     stop('Error: Shape1 parameter should be positive')
  if (min(shape2) <= 0)                     stop('Error: Shape2 parameter should be positive')
  if (min(intensity) < 0)                   stop('Error: Push-intensity parameter should be non-negative')
  if (min(proportion) < 0)                  stop('Error: Push-proportion parameter should be between zero and one')
  if (max(proportion) > 1)                  stop('Error: Push-proportion parameter should be between zero and one')
  
  ######################################################################################################
  
  #Deal with the special case where distribution reduces to the beta distribution
  if (!right) {
    if (intensity == 0) { return(dbeta(x, shape1 = shape1, shape2 = shape2, log = log)) }
    if (proportion == 0)   { return(dbeta(x, shape1 = shape1, shape2 = shape2, log = log)) }
    if (proportion == 1)   { return(dbeta(x, shape1 = shape1, shape2 = shape2 + intensity, log = log)) } }
  if (right) {
    if (intensity == 0) { return(dbeta(1-x, shape1 = shape2, shape2 = shape1, log = log)) }
    if (proportion == 0)   { return(dbeta(1-x, shape1 = shape2, shape2 = shape1, log = log)) }
    if (proportion == 1)   { return(dbeta(1-x, shape1 = shape2 + intensity, shape2 = shape1, log = log)) } }
  
  #Compute log-density vector (left-pushed Gaussian hypergeometric beta distribution)
  if (!right) {
    
    #Compute log-kernel values
    LOGKERNEL <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      xx <- x[i]
      if ((xx == 0)&(shape1 == 1)) { LOGKERNEL[i] <- 0 }
      if ((xx == 1)&(shape2 == 1)) { LOGKERNEL[i] <- intensity*log(1-proportion) }
      if ((xx > 0)&(xx < 1)) { 
        LOGKERNEL[i] <- (shape1-1)*log(xx) + (shape2-1)*log(1-xx) + intensity*log(1-xx*proportion) } }
    
    #Compute scaling constant (using integration)
    LOGSCALE <- scale.pushbeta(upper = 1, shape1 = shape1, shape2 = shape2, intensity = intensity, proportion = proportion, right = right, log = TRUE)
    
    #Generate log-density vector
    LOGDENSITY <- LOGKERNEL - LOGSCALE }
  
  #Compute log-density vector (right-pushed Gaussian hypergeometric beta distribution)
  if (right) {
    
    #Compute log-kernel values
    LOGKERNEL <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      xx <- x[i]
      if ((xx == 0)&(shape1 == 1)) { LOGKERNEL[i] <- intensity*log(1-proportion) }
      if ((xx == 1)&(shape2 == 1)) { LOGKERNEL[i] <- 0 }
      if ((xx > 0)&(xx < 1)) { 
        LOGKERNEL[i] <- (shape1-1)*log(xx) + (shape2-1)*log(1-xx) + intensity*log(1-proportion+xx*proportion) } }
    
    #Compute scaling constant (using integration)
    LOGSCALE <- scale.pushbeta(upper = 1, shape1 = shape1, shape2 = shape2, intensity = intensity, proportion = proportion, right = right, log = TRUE)
    
    #Generate log-density vector
    LOGDENSITY <- LOGKERNEL - LOGSCALE }

  #Return output
  if (log) { LOGDENSITY } else { exp(LOGDENSITY) } }

