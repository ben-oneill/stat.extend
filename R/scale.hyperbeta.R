#' Scaling constant for the Gaussian hypergeometric beta distribution
#'
#' \code{scale.hyperbeta} returns the scaling constant from the distribution.
#'
#' This function computes the scaling constant from the directional form of the Gaussian hypergeometric 
#' beta distribution.  When multiplied by the density kernel it yields the density function.  The scaling 
#' constant is a product of the beta function and Gaussian hypergeometric function, so the present function 
#' is effectively computing the latter over the parameterisation of interest in the distribution.  Further 
#' details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2025) The directional Gaussian hypergeometric beta distributions and their uses in 
#' contaminated binary sampling.
#' 
#' Note: In order to yield greater consistency in integral computation for different upper bounds using the same
#' parameters, the integration method always uses \code{intvals} evenly-spaced points over the unit interval.  If
#' the upper bound of the integral is set to be below one, the method will still use \code{intvals} evenly-spaced
#' points but will ignore the points above the specified upper bound.  This means that if the upper bound of the 
#' integral is set to a low value, you may need to increase the \code{intvals} input in order to use more points 
#' for the integral computation.  The function is vectorised for the \code{upper} input so the user can compute
#' the integral consistently for a set of upper bounds.
#'
#' @usage \code{scale.hyperbeta(upper = 1, shape1, shape2, intensity, proportion, right = FALSE, log = TRUE, warn = TRUE, intvals = 10^6, ...)}
#' @param upper The upper bound of the integral (default is one which gives integral over the entire unit interval)
#' @param shape1 The shape1 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param shape2 The shape2 parameter for the Gaussian hypergeometric beta distribution (positive numeric value)
#' @param intensity The push-intensity parameter for the Gaussian hypergeometric beta distribution (non-negative numeric value)
#' @param proportion The push-proportion parameter for the Gaussian hypergeometric beta distribution (numeric value between zero and one)
#' @param right Logical direction value; \code{TRUE} uses the right-pushed Gaussian hypergeometric beta distribution; \code{FALSE} uses the left-pushed directional Gaussian hypergeometric beta distribution
#' @param log Logical direction value; \code{TRUE} gives the log-integral as the output of the function
#' @param warn Logical direction value; \code{TRUE} means that there will be a warning for low number of points in integral computation
#' @param intvals The number of values over the unit interval to use to approximate the integral
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will be
#' a number giving the scaling constant for the density function (or the logarithm of this constant if preferred)

scale.hyperbeta <- function(upper = 1,
                            shape1, shape2, intensity, proportion, right = FALSE, 
                            log = TRUE, warn = TRUE, intvals = 10^6, ...) {

  #Check that inputs are appropriate type
  if (!is.numeric(upper))                   stop('Error: Input upper is not numeric')
  if (!is.numeric(shape1))                  stop('Error: Shape1 parameter is not numeric')
  if (!is.numeric(shape2))                  stop('Error: Shape2 parameter is not numeric')
  if (!is.numeric(intensity))               stop('Error: Push-intensity parameter is not numeric')
  if (!is.numeric(proportion))              stop('Error: Push-proportion parameter is not numeric')
  if (!is.logical(right))                   stop('Error: Direction input (right) is not a logical value')
  if (!is.logical(log))                     stop('Error: Input log is not a logical value')
  if (!is.numeric(intvals))                 stop('Error: Input intvals is not numeric')
  
  #Check that parameters and options are atomic
  if (length(shape1) != 1)                  stop('Error: Shape1 parameter should be a single number')
  if (length(shape2) != 1)                  stop('Error: Shape2 parameter should be a single number')
  if (length(intensity) != 1)               stop('Error: Push-intensity parameter should be a single number')
  if (length(proportion) != 1)              stop('Error: Push-proportion parameter should be a single number')
  if (length(right) != 1)                   stop('Error: Direction input (right) should be a single logical value')
  if (length(log) != 1)                     stop('Error: Input log should be a single logical value')
  if (length(intvals) != 1)                 stop('Error: Input intvals should be a single number')
  
  #Check that parameters are in allowable range
  if (min(upper) < 0)                       stop('Error: Input upper should be between zero and one')
  if (max(upper) > 1)                       stop('Error: Input upper should be between zero and one')
  if (min(shape1) <= 0)                     stop('Error: Shape1 parameter should be positive')
  if (min(shape2) <= 0)                     stop('Error: Shape2 parameter should be positive')
  if (min(intensity) < 0)                   stop('Error: Push-intensity parameter should be non-negative')
  if (min(proportion) < 0)                  stop('Error: Push-proportion parameter should be between zero and one')
  if (max(proportion) > 1)                  stop('Error: Push-proportion parameter should be between zero and one')
  if (min(intvals) <= 0)                    stop('Error: Input intvals should be greater than zero')
  if (max(intvals) > .Machine$integer.max)  stop('Error: Input intvals is higher than the maximum integer')
  
  ######################################################################################################
  
  #Generate output vector
  n   <- length(upper)
  OUT <- rep(-Inf, n)
  
  #Compute log-scaling constant for the left-pushed Gaussian hypergeometric beta distribution
  if (!right) {
    
    #Set kernel function
    KERNEL <- function(xx) { (xx^(shape1-1))*((1-xx)^(shape2-1))*((1-xx*proportion)^intensity) }
    
    #Compute log-integrals with integrate function
    for (i in 1:n) {
      
      #Compute log-integral using integrate function
      INT <- tryCatch(expr  = {integrate(f = KERNEL, lower = 0, upper = upper[i], rel.tol = .Machine$double.eps)},
                      error = function(e) {NA})
      if (!is.na(INT[1])) {
      if (INT$message == "OK") { OUT[i] <- log(INT$value) } } }
    
    #If integral is computed as zero then this is an error
    #Recompute using importance sampling using evenly-spaced quantiles
    if (any(OUT == -Inf)) {
    
    #Generate quantiles of the beta distribution
    M   <- ceiling(intvals)
    PP  <- (2*(1:M)-1)/(2*M)
    QQ  <- qbeta(PP, shape1 = shape1, shape2 = shape2 + intensity)
    QQQ <- c(0, QQ)
    TT  <- intensity*(log(1-QQ*proportion) - log(1-QQ))
    
    for (i in 1:n) {
    if (OUT[i] == -Inf) {
      
      #Generate log-integral approximation
      WW  <- rep(0, M)
      if (upper[i] < 1) {
      for (k in 1:M) { 
        if (QQQ[k+1] < 1-upper[i]) { WW[k] <- -Inf }
        if ((QQQ[k] <= upper[i])&(QQQ[k+1] > upper[i])) { WW[k] <- log(upper[i]-QQQ[k]) - log(QQQ[k+1]-QQQ[k]) } } }
      
      #Generate output
      OUT[i] <- matrixStats::logSumExp(WW + TT) - log(M) +
                lgamma(shape1) + lgamma(shape2 + intensity) - lgamma(shape1 + shape2 + intensity) } } } }
  
  #Compute moments for the left-pushed Gaussian hypergeometric beta distribution
  if (right) {
    
    #Compute log-integral using integrate function
    KERNEL <- function(xx) { (xx^(shape1-1))*((1-xx)^(shape2-1))*((1-proportion+xx*proportion)^intensity) }
    
    #Compute log-integrals
    for (i in 1:n) {
      
      #Compute log-integral using integrate function
      INT <- tryCatch(expr  = {integrate(f = KERNEL, lower = 0, upper = upper[i], rel.tol = .Machine$double.eps)},
                      error = function(e) {NA})
      if (!is.na(INT[1])) {
      if (INT$message == "OK") { OUT[i] <- log(INT$value) } } }
    
    #If integral is computed as zero then this is an error
    #Recompute using importance sampling using evenly-spaced quantiles
    if (any(OUT == -Inf)) {
      
      #Generate quantiles of the beta distribution
      M   <- ceiling(intvals)
      PP  <- (2*(1:M)-1)/(2*M)
      QQ  <- qbeta(PP, shape1 = shape2, shape2 = shape1 + intensity)
      QQQ <- c(0, QQ)
      TT  <- intensity*(log(1-QQ*proportion) - log(1-QQ))
      
      for (i in 1:n) {
      if (OUT[i] == -Inf) {
        
        #Generate log-integral approximation
        WW <- rep(0, M)
        if (upper[i] < 1) {
          for (k in 1:M) {
            if (QQQ[k+1] < 1-upper[i]) { WW[k] <- -Inf }
            if ((QQQ[k] < 1-upper[i])&(QQQ[k+1] >= 1-upper[i])) { WW[k] <- log(QQQ[k+1]-1+upper[i]) - log(QQQ[k+1]-QQQ[k]) } } }
        OUT[i] <- matrixStats::logSumExp(WW + TT) - log(M) +
                  lgamma(shape1) + lgamma(shape2 + intensity) - lgamma(shape1 + shape2 + intensity) } } } }
  
  #Return the output
  if (log) { OUT } else { exp(OUT) } }
