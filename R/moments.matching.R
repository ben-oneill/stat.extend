#' Moments of the generalised matching distribution
#'
#' \code{moments.match} returns some representative moments from the distribution.
#'
#' This function computes some representative moments from the generalised matching distribution.
#' Further details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2021) A generalised matching distribution for the problem of coincidences.
#'
#' @param size The size parameter for the generalised matching distribution (number of objects to match)
#' @param trials The trials parameter for the generalised matching distribution (number of times the matching game is repeated)
#' @param prob The probability parameter for the generalised matching distribution (probability of known match)
#' @param include.sd Logical value; if \code{TRUE} the output includes the standard deviation
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will be
#' a data frame of moments
#' 
#' @examples 
#' moments.matching(5)

moments.matching <- function(size, trials = 1, prob = 0, include.sd = FALSE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(trials))                  stop('Error: Input trials should be a positive integer')
  if (!is.numeric(size))                    stop('Error: Size parameter is not numeric')
  if (!is.numeric(prob))                    stop('Error: Probability parameter is not numeric')
  if (!is.logical(include.sd))              stop('Error: Input include.sd is not a logical value')

  #Check that parameters are atomic
  if (length(trials) != 1)                  stop('Error: Input trials should be a single positive integer')
  if (length(size)  != 1)                   stop('Error: Size parameter should be a single number')
  if (length(prob)  != 1)                   stop('Error: Probability parameter should be a single number')
  if (length(include.sd)  != 1)             stop('Error: Input include.sd should be a single logical value')

  #Check that parameters are in allowable range
  if (trials == Inf)                        stop('Error: Input trials must be finite')
  m <- as.integer(trials)
  if (m != trials)                          stop('Error: Input trials should be a positive integer')
  if (m <= 0)                               stop('Error: Input trials should be a positive integer')
  if (size < Inf) { n <- as.integer(size) } else { n <- Inf }
  if (n != size)                            stop('Error: Size parameter should be a non-negative integer')
  if (n < 0)                                stop('Error: Size parameter should be a non-negative integer')
  if (prob < 0)                             stop('Error: Probability parameter should be between zero and one')
  if (prob > 1)                             stop('Error: Probability parameter should be between zero and one')

  ########################################################################################
  ################################# TRIVIAL CASES ########################################
  ########################################################################################
  
  #Deal with trivial case where n = 0
  #Distribution is a point-mass on zero
  if (n == 0) {
    MEAN <- 0
    VAR  <- 0
    SKEW <- NA
    KURT <- NA }
  
  #Deal with trivial case where n = 1
  #Distribution is a point-mass on one
  if (n == 1) {
    MEAN <- 1
    VAR  <- 0
    SKEW <- NA
    KURT <- NA }
  
  #Deal with special case where prob = 1
  #Distribution is a point-mass on n
  if (prob == 1) {
    MEAN <- n
    VAR  <- 0
    SKEW <- NA
    KURT <- NA }
  
  #Deal with special case where n = Inf and prob = 0
  #Distribution is the Poisson distribution with unit rate
  if ((n == Inf)&(prob == 0)) {
    MEAN <- 1
    VAR  <- 1
    SKEW <- 1
    KURT <- 4 }
  
  #Deal with special case where n = Inf and prob > 0
  #Distribution is a point-mass on infinity
  if ((n == Inf)&(prob > 0)) {
    MEAN <- Inf
    VAR  <- NA
    SKEW <- NA
    KURT <- NA }
  
  ########################################################################################
  ############################### NON-TRIVIAL CASES ######################################
  ########################################################################################
  
  #Deal with non-trivial case where 1 < n < Inf and prob = 0
  if ((n > 1)&(n < Inf)&(prob == 0)) {
    MEAN <- 1
    VAR  <- 1
    SKEW <- 1*(n>=3)
    KURT <- 1 + 2*(n>=3) + 1*(n>=4) }
  
  #Deal with non-trivial case where 1 < n < Inf and prob > 0
  if ((n > 1)&(n < Inf)&(prob > 0)) {
    T <- prob
    MEAN <- 1 + n*T - T^n
    VAR  <- 1 - T^(2*n) + n*(T-T^2-T^(n-1)-T^n+2*T^(n+1))
    SKEW <- (1 + n*T*(1-3*T+2*T^2) - (1/2)*n*(n-1)*(T^(n-2)) - n*(2*n-1)*(T^(n-1)) + 
             (1/2)*(5*n^2-3*n+2)*(T^n) + 3*n*(n+1)*(T^(n+1)) - 3*n*(n+1)*(T^(n+2)) - 
             3*n*(T^(2*n-1)) - 3*n*(T^(2*n)) + 6*n*(T^(2*n+1)) - 2*(T^(3*n)))/VAR^(3/2)
    KURT <- 3 + (1+n*T-7*n*T^2+12*n*T^3-6*n*T^4 - (1/6)*n*(n-1)*(n-2)*T^(n-3) - 
                 (1/2)*(3*n*(n-1)^2)*T^(n-2) - (1/2)*n*(n^2+3*n-8)*T^(n-1) - 
                 (1/6)*(49*n^3-84*n^2-61*n+6)*T^n - 2*n*(2*n^2-3*n)*T^(n+1) - 
                 6*n*(n^2+3*n+2)*T^(n+2) + 4*n*(n+1)*(n+2)*T^(n+3) - n*(5*n-2)*T^(2*n-2) -
                 2*n*(7*n-2)*T^(2*n-1) - (n^2-12*n-2)*T^(2*n) + 12*n*(2*n+1)*T^(2*n+1) - 
                 12*n*(2*n+1)*T^(2*n+2) - 12*n*T^(3*n) + 24*n*T^(3*n+1) - 6*T^(4*n))/VAR^2 }
  
  #Determine moments for multiple trials
  if (m > 1) {
    MEAN.TOTAL <- m*MEAN
     VAR.TOTAL <- m*VAR
    SKEW.TOTAL <- SKEW/sqrt(m)
    KURT.TOTAL <- 3 + (KURT-3)/m }
  else {
    MEAN.TOTAL <- MEAN
    VAR.TOTAL <- VAR
    SKEW.TOTAL <- SKEW
    KURT.TOTAL <- KURT}

  #Generate output
  if (include.sd) {
    OUT <- data.frame(mean = 0, var = 0, sd = 0, skew = 0, kurt = 0, excess.kurt = 0) } else {
    OUT <- data.frame(mean = 0, var = 0, skew = 0, kurt = 0, excess.kurt = 0) }
  OUT$mean <- MEAN.TOTAL
  OUT$var  <- VAR.TOTAL
  if (include.sd) { OUT$sd <- sqrt(VAR.TOTAL) }
  OUT$skew <- SKEW.TOTAL
  OUT$kurt <- KURT.TOTAL
  OUT$excess.kurt <- KURT.TOTAL - 3

  #Return the output
  OUT }
