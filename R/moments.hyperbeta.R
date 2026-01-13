#' Moments of the Gaussian hypergeometric beta distribution
#'
#' \code{moments.hyperbeta} returns some representative moments from the distribution.
#'
#' This function computes some representative moments from the directional form of the Gaussian 
#' hypergeometric beta distribution.  Further details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2025) Directional Gaussian hypergeometric beta distributions and their uses in
#' contaminated binary sampling.
#'
#' @usage \code{moments.hyperbeta(shape1, shape2, intensity, proportion, right = FALSE, include.sd = FALSE)}
#' @param shape1 The shape1 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param shape2 The shape2 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param intensity The push-intensity parameter for the Gaussian hypergeometric beta distribution (non-negative numeric value)
#' @param proportion The push-proportion parameter for the Gaussian hypergeometric beta distribution (numeric value between zero and one)
#' @param right Logical direction value; \code{TRUE} uses the right-pushed Gaussian hypergeometric beta distribution; \code{FALSE} uses the left-pushed Gaussian hypergeometric beta distribution
#' @param include.sd Logical value; if \code{TRUE} the output includes the standard deviation
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will be
#' a data frame of moments

moments.hyperbeta <- function(shape1, shape2, intensity, proportion, right = FALSE, include.sd = FALSE) {

  #Check that inputs are appropriate type
  if (!is.numeric(shape1))                  stop('Error: Shape1 parameter is not numeric')
  if (!is.numeric(shape2))                  stop('Error: Shape2 parameter is not numeric')
  if (!is.numeric(intensity))               stop('Error: Push-intensity parameter is not numeric')
  if (!is.numeric(proportion))              stop('Error: Push-proportion parameter is not numeric')
  if (!is.logical(right))                   stop('Error: Direction input (right) is not a logical value')
  if (!is.logical(include.sd))              stop('Error: include.sd option is not a logical value')
  
  #Check that parameters and options are atomic
  if (length(shape1) != 1)                  stop('Error: Shape1 parameter should be a single number')
  if (length(shape2) != 1)                  stop('Error: Shape2 parameter should be a single number')
  if (length(intensity) != 1)               stop('Error: Push-intensity parameter should be a single number')
  if (length(proportion) != 1)              stop('Error: Push-proportion parameter should be a single number')
  if (length(right) != 1)                   stop('Error: Direction input (right) should be a single logical value')
  if (length(include.sd) != 1)              stop('Error: include.sd option should be a single logical value')
  
  #Check that parameters are in allowable range
  if (min(shape1) <= 0)                     stop('Error: Shape1 parameter should be positive')
  if (min(shape2) <= 0)                     stop('Error: Shape2 parameter should be positive')
  if (min(intensity) < 0)                   stop('Error: Push-intensity parameter should be non-negative')
  if (min(proportion) < 0)                  stop('Error: Push-proportion parameter should be between zero and one')
  if (max(proportion) > 1)                  stop('Error: Push-proportion parameter should be between zero and one')
  
  ######################################################################################################
  
  #Compute the raw moments
  LOGSCALE <- scale.hyperbeta(upper = 1, shape1 = shape1,   shape2 = shape2, intensity = intensity, proportion = proportion, right, log = TRUE)
  LOGINT1  <- scale.hyperbeta(upper = 1, shape1 = shape1+1, shape2 = shape2, intensity = intensity, proportion = proportion, right, log = TRUE)
  LOGINT2  <- scale.hyperbeta(upper = 1, shape1 = shape1+2, shape2 = shape2, intensity = intensity, proportion = proportion, right, log = TRUE)
  LOGINT3  <- scale.hyperbeta(upper = 1, shape1 = shape1+3, shape2 = shape2, intensity = intensity, proportion = proportion, right, log = TRUE)
  LOGINT4  <- scale.hyperbeta(upper = 1, shape1 = shape1+4, shape2 = shape2, intensity = intensity, proportion = proportion, right, log = TRUE)
  RAW1     <- exp(LOGINT1 - LOGSCALE)
  RAW2     <- exp(LOGINT2 - LOGSCALE)
  RAW3     <- exp(LOGINT3 - LOGSCALE)
  RAW4     <- exp(LOGINT4 - LOGSCALE)
  
  #Compute the central moments
  CENT2 <- RAW2 - RAW1^2
  CENT3 <- RAW3 - 3*RAW1*RAW2 + 2*RAW1^3
  CENT4 <- RAW4 - 4*RAW1*RAW3 + 6*(RAW1^2)*RAW2 - 3*RAW1^4
  
  #Compute the moments
  MEAN <- RAW1
  VAR  <- CENT2
  SKEW <- CENT3/(CENT2^(3/2))
  KURT <- CENT4/(CENT2^2)
  
  #Generate output
  if (include.sd) {
    OUT <- data.frame(mean = MEAN, var = VAR, sd = sqrt(VAR), skew = SKEW, kurt = KURT, excess.kurt = KURT - 3) 
  } else {
    OUT <- data.frame(mean = MEAN, var = VAR, skew = SKEW, kurt = KURT, excess.kurt = KURT - 3) }
  
  #Return the output
  OUT }
