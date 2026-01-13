#' Quantile function of the Gaussian hypergeometric beta distribution
#'
#' \code{qhyperbeta} returns the quantile values for the arguments.
#'
#' This function computes quantiles from the quantile function of the directional form of the Gaussian
#' hypergeometric beta distribution.  Further details on the distribution can be found in the following 
#' paper:
#'
#' O'Neill, B. (2025) Directional Gaussian hypergeometric beta distributions and their uses in 
#' contaminated binary sampling.
#'
#' @usage \code{qhyperbeta(p, shape1, shape2, intensity, proportion, right = FALSE, lower.tail = TRUE, log.p = FALSE)}
#' @param p A vector of numeric values to be used as probability/log-probability arguments for the quantile function
#' @param shape1 The shape1 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param shape2 The shape2 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param intensity The push-shape parameter for the Gaussian hypergeometric beta distribution (non-negative numeric value)
#' @param proportion The push-proportion parameter for the Gaussian hypergeometric beta distribution (numeric value between zero and one)
#' @param right Logical direction value; \code{TRUE} uses the right-pushed Gaussian hypergeometric beta distribution; \code{FALSE} uses the left-pushed Gaussian hypergeometric beta distribution
#' @param lower.tail A logical value specifying whether results are from the cumulative distribution function
#' or the corresponding survival function
#' @param log.p A logical value specifying whether input arguments are log-probabilities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range and arguments are numeric)
#' then the output will be a vector of quantiles corresponding to the vector argument p

qhyperbeta <- function(p, shape1, shape2, intensity, proportion, right = FALSE, lower.tail = TRUE, log.p = FALSE) {

  #Check that inputs are appropriate type
  if (!is.numeric(p))                       stop('Error: Argument p is not numeric')
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
  if (length(proportion)!= 1)               stop('Error: Push-proportion parameter should be a single number')
  if (length(right) != 1)                   stop('Error: Direction input (right) should be a single logical value')
  if (length(lower.tail) != 1)              stop('Error: lower.tail option should be a single logical value')
  if (length(log.p) != 1)                   stop('Error: log.p option should be a single logical value')
  
  #Check the input p
  if (log.p) {
    if (max(p) > 0)                         stop('Error: Input log-probabilities in p cannot be positive')
    LOG.P <- p
  } else {
    if (min(p) < 0)                         stop('Error: Input probabilities in p cannot be negative')
    if (max(p) > 1)                         stop('Error: Input probabilities in p cannot be greater than one')
    LOG.P <- log(p) }
  
  #Check that parameters are in allowable range
  if (min(shape1) <= 0)                     stop('Error: Shape1 parameter should be positive')
  if (min(shape2) <= 0)                     stop('Error: Shape2 parameter should be positive')
  if (min(intensity) < 0)                   stop('Error: Push-intensity parameter should be non-negative')
  if (min(proportion) < 0)                  stop('Error: Push-proportion parameter should be between zero and one')
  if (max(proportion) > 1)                  stop('Error: Push-proportion parameter should be between zero and one')
  
  ######################################################################################################
  
  #Deal with the special case where distribution reduces to the beta distribution
  if (!right) {
    if (intensity == 0)    { return(qbeta(LOG.P, shape1 = shape1, shape2 = shape2, lower.tail = lower.tail, log.p = TRUE)) }
    if (proportion == 0)   { return(qbeta(LOG.P, shape1 = shape1, shape2 = shape2, lower.tail = lower.tail, log.p = TRUE)) }
    if (proportion == 1)   { return(qbeta(LOG.P, shape1 = shape1, shape2 = shape2 + intensity, lower.tail = lower.tail, log.p = TRUE)) } }
  if (right) {
    LOG.P <- VGAM::log1mexp(-LOG.P)
    if (intensity == 0)    { return(qbeta(LOG.P, shape1 = shape2, shape2 = shape1, lower.tail = !lower.tail, log.p = TRUE)) }
    if (proportion == 0)   { return(qbeta(LOG.P, shape1 = shape2, shape2 = shape1, lower.tail = !lower.tail, log.p = TRUE)) }
    if (proportion == 1)   { return(qbeta(LOG.P, shape1 = shape2 + intensity, shape2 = shape1, lower.tail = !lower.tail, log.p = TRUE)) } }
  
  #Compute quantile vector (left-pushed Gaussian hypergeometric beta distribution)
  if (!right) {
    
    #Compute scaling constant (using integration)
    LOGSCALE <- scale.hyperbeta(upper = 1, shape1 = shape1, shape2 = shape2, intensity = intensity, proportion = proportion, right = right, log = TRUE)
    
    #Compute quantiles
    QUANTILES <- rep(0, length(LOG.P))
    if (!lower.tail) { LOG.P <- VGAM::log1mexp(-LOG.P) }
    for (i in 1:length(LOG.P)) {
      
      if (LOG.P[i] == -Inf) { QUANTILES[i] <- 0 }
      if (LOG.P[i] == 0)    { QUANTILES[i] <- 1 }
      if ((LOG.P[i] > -Inf)&(LOG.P[i] < 0)) {
        
        #Compute quantiles from nonlinear optimisation
        kk <- LOG.P[i]
        q0 <- qbeta(LOG.P[i], shape1 = shape1, shape2 = shape2 + proportion*intensity, lower.tail = TRUE, log.p = TRUE)
        x0 <- log(q0/(1-q0))
        OBJECTIVE <- function(xx) {
          qq  <- exp(xx)/(1+exp(xx))
          LOGPROB <- scale.hyperbeta(upper = qq, shape1 = shape1, shape2 = shape2, intensity = intensity, proportion = proportion, right = right, log = TRUE) - LOGSCALE
          (LOGPROB - kk)^2 }
        NLM <- nlm(f = OBJECTIVE, p = x0, gradtol = 10e-21)
        XX <- NLM$estimate
        QUANTILES[i] <- exp(XX)/(1+exp(XX)) } } }
  
  #Compute quantile vector (right-pushed Gaussian hypergeometric beta distribution)
  if (right) {
    
    #Compute scaling constant (using integration)
    LOGSCALE <- scale.hyperbeta(upper = 1, shape1 = shape1, shape2 = shape2, intensity = intensity, proportion = proportion, right = right, log = TRUE)
    
    #Compute quantiles
    QUANTILES <- rep(0, length(LOG.P))
    if (!lower.tail) { LOG.P <- VGAM::log1mexp(-LOG.P) }
    for (i in 1:length(LOG.P)) {
      
      if (LOG.P[i] == -Inf) { QUANTILES[i] <- 0 }
      if (LOG.P[i] == 0)    { QUANTILES[i] <- 1 }
      if ((LOG.P[i] > -Inf)&(LOG.P[i] < 0)) {
        
        #Compute quantiles from nonlinear optimisation
        kk <- LOG.P[i]
        q0 <- qbeta(LOG.P[i], shape1 = shape1, shape2 = shape2 + proportion*intensity, lower.tail = TRUE, log.p = TRUE)
        x0 <- log(q0/(1-q0))
        OBJECTIVE <- function(xx) {
          qq  <- exp(xx)/(1+exp(xx))
          LOGPROB <- scale.hyperbeta(upper = qq, shape1 = shape1, shape2 = shape2, intensity = intensity, proportion = proportion, right = right, log = TRUE) - LOGSCALE
          (LOGPROB - kk)^2 }
        NLM <- nlm(f = OBJECTIVE, p = x0, gradtol = 10e-21)
        XX <- NLM$estimate
        QUANTILES[i] <- exp(XX)/(1+exp(XX)) } } }

  #Return output
  QUANTILES }
