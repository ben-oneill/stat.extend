#' Random generation from the Gaussian hypergeometric beta distribution
#'
#' \code{rhyperbeta} returns random variables from the distribution.
#'
#' This function generates random variables from the directional form of the Gaussian hypergeometric 
#' beta distribution.  Further details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2025) Directional Gaussian hypergeometric beta distributions and their uses in
#' contaminated binary sampling.
#'
#' @usage \code{rhyperbeta(n, shape1, shape2, intensity, proportion, right = FALSE)}
#' @param n The number of random variables to generate
#' @param shape1 The shape1 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param shape2 The shape2 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param intensity The push-intensity parameter for the Gaussian hypergeometric beta distribution (non-negative numeric value)
#' @param proportion The push-proportion parameter for the Gaussian hypergeometric beta distribution (numeric value between zero and one)
#' @param right Logical direction value; \code{TRUE} uses the right-pushed Gaussian hypergeometric beta distribution; \code{FALSE} uses the left-pushed Gaussian hypergeometric beta distribution
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range and arguments are numeric)
#' then the output will be a vector of random variables from the distribution

rhyperbeta <- function(n, shape1, shape2, intensity, proportion, right = FALSE) {

  #Check that inputs are appropriate type
  if (!is.numeric(n))                       stop('Error: Argument n is not numeric')
  if (!is.numeric(shape1))                  stop('Error: Shape1 parameter is not numeric')
  if (!is.numeric(shape2))                  stop('Error: Shape2 parameter is not numeric')
  if (!is.numeric(intensity))               stop('Error: Push-intensity parameter is not numeric')
  if (!is.numeric(proportion))              stop('Error: Push-proportion parameter is not numeric')
  if (!is.logical(right))                   stop('Error: Direction input (right) is not a logical value')
  
  #Check that parameters and options are atomic
  if (length(n) != 1)                       stop('Error: Argument n should be a single number')
  if (length(shape1) != 1)                  stop('Error: Shape1 parameter should be a single number')
  if (length(shape2) != 1)                  stop('Error: Shape2 parameter should be a single number')
  if (length(intensity) != 1)               stop('Error: Push-intensity parameter should be a single number')
  if (length(proportion) != 1)              stop('Error: Push-proportion parameter should be a single number')
  if (length(right) != 1)                   stop('Error: Direction input (right) should be a single logical value')
  
  #Check input n
  nn <- as.integer(n)
  if (nn != n)                              stop('Error: Argument n should be an integer')
  if (nn < 0)                               stop('Error: Argument n should be a non-negative integer')
  if (nn == Inf)                            stop('Error: Argument n should be finite')
  
  #Check that parameters are in allowable range
  if (min(shape1) <= 0)                     stop('Error: Shape1 parameter should be positive')
  if (min(shape2) <= 0)                     stop('Error: Shape2 parameter should be positive')
  if (min(intensity) < 0)                   stop('Error: Push-intensity parameter should be non-negative')
  if (min(proportion) < 0)                  stop('Error: Push-proportion parameter should be between zero and one')
  if (max(proportion) > 1)                  stop('Error: Push-proportion parameter should be between zero and one')
  
  ######################################################################################################
  
  #Deal with the special case where distribution reduces to the beta distribution
  if (!right) {
    if (intensity == 0)  { return(rbeta(n, shape1 = shape1, shape2 = shape2)) }
    if (proportion == 0) { return(rbeta(n, shape1 = shape1, shape2 = shape2)) }
    if (proportion == 1) { return(rbeta(n, shape1 = shape1, shape2 = shape2 + intensity)) } }
  if (right) {
    if (intensity == 0)  { return(1-rbeta(n, shape1 = shape2, shape2 = shape1)) }
    if (proportion == 0) { return(1-rbeta(n, shape1 = shape2, shape2 = shape1)) }
    if (proportion == 1) { return(1-rbeta(n, shape1 = shape2 + intensity, shape2 = shape1)) } }
  
  #Generate random variables (left-pushed Gaussian hypergeometric beta distribution)
  if (!right) {
    
    #Compute scaling constant (using integration)
    LOGSCALE <- scale.hyperbeta(upper = 1, shape1 = shape1, shape2 = shape2, intensity = intensity, proportion = proportion, right = right, log = TRUE)
    
    #Compute quantiles
    OUT   <- rep(0, n)
    LOG.P <- log(runif(n))
    for (i in 1:n) {
      
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
      OUT[i] <- exp(XX)/(1+exp(XX)) } }

  #Compute quantile vector (right-pushed Gaussian hypergeometric beta distribution)
  if (right) {
    
    #Compute scaling constant (using integration)
    LOGSCALE <- scale.hyperbeta(upper = 1, shape1 = shape1, shape2 = shape2, intensity = intensity, proportion = proportion, right = right, log = TRUE)
    
    #Compute quantiles
    OUT   <- rep(0, n)
    LOG.P <- log(runif(n))
    for (i in 1:n) {
      
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
      OUT[i] <- exp(XX)/(1+exp(XX)) } }

  #Return output
  OUT }
