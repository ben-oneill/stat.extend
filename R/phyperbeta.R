#' Cumulative distribution function of the Gaussian hypergeometric beta distribution
#'
#' \code{phyperbeta} returns the probability or log-probability values for the arguments.
#'
#' This function computes probabilities or log-probabilities from the cumulative distribution 
#' function of the directional form of the Gaussian hypergeometric beta distribution.  Further 
#' details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2025) Directional Gaussian hypergeometric beta distributions and their uses in 
#' contaminated binary sampling.
#'
#' @usage \code{phyperbeta(x, shape1, shape2, intensity, proportion, right = FALSE, lower.tail = TRUE, log.p = FALSE)}
#' @param x A vector of numeric values to be used as arguments for the cumulative distribution function
#' @param shape1 The shape1 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param shape2 The shape2 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param intensity The push-intensity parameter for the Gaussian hypergeometric beta distribution (non-negative numeric value)
#' @param proportion The push-proportion parameter for the Gaussian hypergeometric beta distribution (numeric value between zero and one)
#' @param right Logical direction value; \code{TRUE} uses the right-pushed Gaussian hypergeometric beta distribution; \code{FALSE} uses the left-pushed Gaussian hypergeometric beta distribution
#' @param lower.tail A logical value specifying whether results are from the cumulative distribution function
#' or the corresponding survival function
#' @param log.p A logical value specifying whether results should be returned as log-probabilities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range and arguments are numeric)
#' then the output will be a vector of probabilities/log-probabilities corresponding to the vector argument x

phyperbeta <- function(x, shape1, shape2, intensity, proportion, right = FALSE, lower.tail = TRUE, log.p = FALSE) {

  #Check that inputs are appropriate type
  if (!is.numeric(x))                       stop('Error: Argument x is not numeric')
  if (!is.numeric(shape1))                  stop('Error: Shape1 parameter is not numeric')
  if (!is.numeric(shape2))                  stop('Error: Shape2 parameter is not numeric')
  if (!is.numeric(intensity))               stop('Error: Push-intensity parameter is not numeric')
  if (!is.numeric(proportion))              stop('Error: Push-proportion parameter is not numeric')
  if (!is.logical(right))                   stop('Error: Direction input (right) is not a logical value')
  if (!is.logical(lower.tail))              stop('Error: lower.tail option is not a logical value')
  if (!is.logical(log.p))                   stop('Error: log.p option is not a logical value')
  
  #Check that parameters and options are atomic
  if (length(shape1) != 1)                  stop('Error: Shape1 parameter should be a single number')
  if (length(shape2) != 1)                  stop('Error: Shape2 parameter should be a single number')
  if (length(intensity) != 1)               stop('Error: Push-intensity parameter should be a single number')
  if (length(proportion) != 1)              stop('Error: Push-proportion parameter should be a single number')
  if (length(right) != 1)                   stop('Error: Direction input (right) should be a single logical value')
  if (length(lower.tail) != 1)              stop('Error: lower.tail option should be a single logical value')
  if (length(log.p) != 1)                   stop('Error: log.p option should be a single logical value')
  
  #Check that parameters are in allowable range
  if (min(shape1) <= 0)                     stop('Error: Shape1 parameter should be positive')
  if (min(shape2) <= 0)                     stop('Error: Shape2 parameter should be positive')
  if (min(intensity) < 0)                   stop('Error: Push-intensity parameter should be non-negative')
  if (min(proportion) < 0)                  stop('Error: Push-proportion parameter should be between zero and one')
  if (max(proportion) > 1)                  stop('Error: Push-proportion parameter should be between zero and one')
  
  ######################################################################################################
  
  #Deal with the special case where distribution reduces to the beta distribution
  if (!right) {
    if (intensity == 0)    { return(pbeta(x, shape1 = shape1, shape2 = shape2, lower.tail = lower.tail, log.p = log.p)) }
    if (proportion == 0)   { return(pbeta(x, shape1 = shape1, shape2 = shape2, lower.tail = lower.tail, log.p = log.p)) }
    if (proportion == 1)   { return(pbeta(x, shape1 = shape1, shape2 = shape2 + intensity, lower.tail = lower.tail, log.p = log.p)) } }
  if (right) {
    if (intensity == 0)    { return(pbeta(1-x, shape1 = shape2, shape2 = shape1, lower.tail = !lower.tail, log.p = log.p)) }
    if (proportion == 0)   { return(pbeta(1-x, shape1 = shape2, shape2 = shape1, lower.tail = !lower.tail, log.p = log.p)) }
    if (proportion == 1)   { return(pbeta(1-x, shape1 = shape2 + intensity, shape2 = shape1, lower.tail = !lower.tail, log.p = log.p)) } }
  
  #Compute log-probability vector
  LOGSCALE <- scale.hyperbeta(upper = 1, shape1 = shape1, shape2 = shape2, intensity = intensity, proportion = proportion, right = right, log = TRUE)
  LOGPROB  <- scale.hyperbeta(upper = pmax(0, pmin(x, 1)), shape1 = shape1, shape2 = shape2, intensity = intensity, proportion = proportion, right = right, log = TRUE) - LOGSCALE
  
  #Flip probabilities if using upper tail
  if (!lower.tail) { LOGPROB <- VGAM::log1mexp(-LOGPROB) }

  #Return output
  if (log.p) { LOGPROB } else { exp(LOGPROB) } }
