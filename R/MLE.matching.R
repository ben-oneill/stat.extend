#' Maximum likelihood estimator (MLE) in the generalised matching distribution
#'
#' \code{MLE.matching} returns the maximum likelihood estimator (MLE) for the data.
#'
#' This function computes the maximum likelihood estimator (MLE) from data consisting of IID samples
#' from the generalised matching distribution.  Further details on the distribution can be found in 
#' the following paper:
#'
#' @references 
#' O'Neill, B. (2021) A generalised matching distribution for the problem of coincidences.
#'
#' @param x A vector of numeric values to be used as arguments for the mass function
#' @param size The size parameter for the generalised matching distribution (number of objects to match)
#' @param CI.method The method used to compute the confidence interval ('asymptotic' or 'bootstrap')
#' @param conf.level The width of the CI
#' @param bootstrap.sims The number of bootstrap simulations used in the bootstrap confidence interval
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will
#' be a list of outputs for the MLE
#' 
#' @examples 
#' 
#' X <- rmatching(20, 5, prob=.1)
#' 
#' # For comparison
#' # MASS::fitdistr(X, dmatching, start=list(prob=.5), size=5, lower=c(prob=0), upper=c(prob=1))
#' 
#' MLE.matching(X, 5)
#' 

MLE.matching <- function(x, size, CI.method = 'asymptotic', conf.level = 0.95, bootstrap.sims = 10^3) {

  #Extract data name
  DATA.NAME <- deparse(substitute(x))
  
  #Check inputs x and size
  if (!is.numeric(x))                       stop('Error: Input x should be numeric')
  if (!is.numeric(size))                    stop('Error: Size parameter should be a non-negative integer')
  if (length(size) != 1)                    stop('Error: Size parameter should be a single non-negative integer')
  if (size < Inf) { n <- as.integer(size) } else { n <- Inf }
  if (n != size)                            stop('Error: Size parameter should be a non-negative integer')
  if (n < 0)                                stop('Error: Size parameter should be a non-negative integer')
  m <- length(x)
  for (i in 1:m) {
    if (!(x[i] %in% 0:n)) stop(paste0('Element ', i, ' of the data vector x is not in the support'))
    if (x[i] == n-1)      stop(paste0('Element ', i, ' of the data vector x is not in the support')) }
  
  #Check other inputs
  if (!(CI.method %in% c('asymptotic', 'bootstrap')))   stop('Error: CI.method not recognised')
  if (!is.numeric(conf.level))              stop('Error: Input conf.level must be numeric')
  if (length(conf.level) != 1)              stop('Error: Input conf.level must be a single numeric value')
  if (min(conf.level) < 0)                  stop('Error: Input conf.level cannot be less than zero')
  if (max(conf.level) > 1)                  stop('Error: Input conf.level cannot be greater than one')
  if (!is.numeric(bootstrap.sims))          stop('Error: Input bootstrap.sims must be a positive integer')
  if (length(bootstrap.sims) != 1)          stop('Error: Input bootstrap.sims must be a single positive integer')
  SIMS <- as.integer(bootstrap.sims)
  if (bootstrap.sims != SIMS)               stop('Error: Input bootstrap.sims must be a positive integer')
  if (min(bootstrap.sims) < 1)              stop('Error: Input bootstrap.sims must be a positive integer')
  
  ########################################################################################
  ################################# TRIVIAL CASES ########################################
  ########################################################################################

  #Deal with trivial case where n = 0 or n = 1
  #Distribution is a point-mass distribution (MLE is NA)
  if ((n == 0)|(n == 1)) {
    return(list(MLE.prob = NA, standard.error = NA,
                maxloglike = 0, maxloglike.mean = 0, CI.prob = NA)) }

  ########################################################################################
  ############################### NON-TRIVIAL CASES ######################################
  ########################################################################################

  #Compute a matrix of classical log-probability values
  BASE.LOGPROBS <- matrix(-Inf, nrow = n+1, ncol = n+1)
  rownames(BASE.LOGPROBS) <- sprintf('size[%s]',  0:n)
  colnames(BASE.LOGPROBS) <- sprintf('match[%s]', 0:n)
  BASE.LOGPROBS[1,1] <- 0
  for (nn in 1:n) {
    BASE.LOGPROBS[nn+1, nn+1] <- -lfactorial(nn)
    for (i in 1:nn) {
      k <- nn-i
      T1 <- log(nn-k-1) + BASE.LOGPROBS[nn+1, k+2]
      T2 <- ifelse(k < nn-1, log(k+2) + BASE.LOGPROBS[nn+1, k+3], -Inf)
      BASE.LOGPROBS[nn+1, k+1] <- log(k+1) - log(nn-k) + matrixStats::logSumExp(c(T1, T2)) }
    BASE.LOGPROBS[nn+1, ] <- BASE.LOGPROBS[nn+1, ] - matrixStats::logSumExp(BASE.LOGPROBS[nn+1, ]) }
  
  #Set matrix M
  M <- matrix(-Inf, nrow = n+1, ncol = n+1)
  rownames(M) <- sprintf('k[%s]', 0:n)
  colnames(M) <- sprintf('l[%s]', 0:n)
  for (k in 0:n) {
    for (l in 0:k) {
      M[k+1, l+1] <- BASE.LOGPROBS[n-l+1, k-l+1] } }

  #Set objective function
  NEGLOGLIKE <- function(phi, hessian.prob = FALSE) {
    
    #Set probability parameter
    prob <- exp(phi)/(exp(phi) + exp(-phi))
    
    #Compute vector of log-probability values
    LOGPROBS <- rep(-Inf, n+1)
    LOGBINOM <- dbinom(0:n, size = n, prob = prob, log = TRUE)
    GRAD  <- rep(0, n+1)
    HESS  <- rep(0, n+1)
    T1 <- 2*((0:k) - n*prob)
    T2 <- 8*(1/2 - prob)*((0:k) - n*prob)  + 4*((0:k-1)*(0:k) - 2*(0:k)*(n-1)*prob + n*(n-1)*prob^2)
    for (k in 0:n) {
      LOGPROBS[k+1] <- matrixStats::logSumExp(LOGBINOM + M[k+1, ]) }
    LOGPROBS <- LOGPROBS - matrixStats::logSumExp(LOGPROBS)
    for (k in 0:n) {
      NUM1 <- sum(T1*exp(LOGBINOM + M[k+1, ]))
      NUM2 <- sum(T2*exp(LOGBINOM + M[k+1, ]))
      GRAD[k+1]  <- NUM1/exp(LOGPROBS[k+1])
      HESS[k+1]  <- NUM2/exp(LOGPROBS[k+1]) -  GRAD[k+1]^2 }
    
    #Compute log-likelihood
    LLL <- rep(-Inf, m)
    GGG <- rep(-Inf, m)
    HHH <- rep(-Inf, m)
    for (i in 1:m) {
      xx <- x[i]
      LLL[i] <- - LOGPROBS[xx+1]
      GGG[i] <- - GRAD[xx+1]
      HHH[i] <- - HESS[xx+1] }
    OUT <- sum(LLL)
    attr(OUT, 'gradient') <- sum(GGG)
    attr(OUT, 'hessian')  <- sum(HHH)
    
    #Return output
    OUT }
  
  #Set starting point for iteration
  start.prob <- (mean(x)+1)/(n+1)
  start.phi  <- atanh(2*start.prob - 1)
  
  #Compute MLE
  NLM <- nlm(f = NEGLOGLIKE, p = start.phi, hessian = TRUE, gradtol = 1e-20, steptol = 1e-20)
  MLE.phi  <- NLM$estimate
  MLE.prob <- exp(MLE.phi)/(exp(MLE.phi) + exp(-MLE.phi))
  SE.phi   <- 1/sqrt(c(NLM$hessian))
  SE.prob  <- 2*MLE.prob*(1-MLE.prob)/sqrt(NLM$hessian)
  MLE.DF   <- data.frame(prob = c(MLE.prob, SE.prob), phi = c(MLE.phi, SE.phi))
  rownames(MLE.DF) <- c('MLE', 'SE.asymptotic')
  
  #Compute maximum log-likelihood
  LL.MAX   <- - NLM$minimum
  
  #Generate confidence interval
  alpha     <- 1-conf.level
  alpha.low <- (max(mean(x)-1, 0)/(n-1))*alpha
  LOW.phi   <- qnorm(alpha.low,         mean = MLE.phi, sd = SE.phi)
  HIGH.phi  <- qnorm(1-alpha+alpha.low, mean = MLE.phi, sd = SE.phi)
  LOW.prob  <- (1 + tanh(LOW.phi))/2
  HIGH.prob <- (1 + tanh(HIGH.phi))/2
  CI.PROB   <- c(LOW.prob, HIGH.prob)
  attr(CI.PROB, 'conf.level') <- conf.level
  attr(CI.PROB, 'method')     <- CI.method
  
  #Set output
  OUT.MLE <- list(size = n, data.name = DATA.NAME, data = x,
                  MLE = MLE.DF, maxloglike = LL.MAX, maxloglike.mean = LL.MAX/m,
                  code = NLM$code, iterations = NLM$iterations, CI.prob = CI.PROB)
  
  #Create bootstrap CI (if selected)
  if (CI.method == 'bootstrap') {
    
    #Generate bootstrap MLEs
    MLE.phi.bootstrap <- rep(0, bootstrap.sims)
    for (j in 1:bootstrap.sims) {
      
      #Create bootstrap sample
      SAMPLE <- sample(x, size = m, replace = TRUE)
      
      #Set objective function
      NEGLOGLIKE.bs <- function(phi) {
        
        #Set probability parameter
        prob <- exp(phi)/(exp(phi) + exp(-phi))
        
        #Compute vector of log-probability values
        LOGPROBS <- rep(-Inf, n+1)
        LOGBINOM <- dbinom(0:n, size = n, prob = prob, log = TRUE)
        GRAD <- rep(0, n+1)
        HESS <- rep(0, n+1)
        T1 <- 2*((0:k) - n*prob)
        T2 <- 8*(1/2 - prob)*((0:k) - n*prob)  + 4*((0:k-1)*(0:k) - 2*(0:k)*(n-1)*prob + n*(n-1)*prob^2)
        for (k in 0:n) {
          LOGPROBS[k+1] <- matrixStats::logSumExp(LOGBINOM + M[k+1, ]) }
        LOGPROBS <- LOGPROBS - matrixStats::logSumExp(LOGPROBS)
        for (k in 0:n) {
          NUM1 <- sum(T1*exp(LOGBINOM + M[k+1, ]))
          NUM2 <- sum(T2*exp(LOGBINOM + M[k+1, ]))
          GRAD[k+1] <- NUM1/exp(LOGPROBS[k+1])
          HESS[k+1] <- NUM2/exp(LOGPROBS[k+1]) - GRAD[k+1]^2 }
        
        #Compute log-likelihood
        LLL <- rep(-Inf, m)
        GGG <- rep(-Inf, m)
        HHH <- rep(-Inf, m)
        for (i in 1:m) {
          xx <- SAMPLE[i]
          LLL[i] <- - LOGPROBS[xx+1]
          GGG[i] <- - GRAD[xx+1]
          HHH[i] <- - HESS[xx+1] }
        OUT <- sum(LLL)
        attr(OUT, 'gradient') <- sum(GGG)
        attr(OUT, 'hessian')  <- sum(HHH)
     
        #Return output
        OUT }
      
      #Set starting point for iteration
      start.prob <- (mean(SAMPLE)+1)/(n+1)
      start.phi  <- atanh(2*start.prob - 1)
      
      #Compute bootstrap MLE
      NLM <- nlm(f = NEGLOGLIKE.bs, p = start.phi, hessian = TRUE, gradtol = 1e-20, steptol = 1e-20)
      MLE.phi.bootstrap[j] <- NLM$estimate }
    
    #Update the MLE in output
    MLE.phi.bootstrap  <- sort(MLE.phi.bootstrap)
    MLE.prob.bootstrap <- exp(MLE.phi.bootstrap)/(exp(MLE.phi.bootstrap) + exp(-MLE.phi.bootstrap))
    SE.phi  <- sd(MLE.phi.bootstrap)
    SE.prob <- sd(MLE.prob.bootstrap)
    MLE.DF[2,] <- c(SE.prob, SE.phi)
    rownames(MLE.DF) <- c('MLE', 'SE.bootstrap')
    OUT.MLE$MLE <- MLE.DF
    
    #Generate confidence interval
    VALUES    <- c(0, rep(MLE.prob.bootstrap, each = 2), 1)
    alpha     <- 1-conf.level
    alpha.low <- (max(mean(x)-1, 0)/(n-1))*alpha
    LOW.prob  <- quantile(x = VALUES, prob = alpha.low)
    HIGH.prob <- quantile(x = VALUES, prob = 1-alpha+alpha.low)
    CI.PROB   <- c(LOW.prob, HIGH.prob)
    names(CI.PROB) <- c('Low', 'High')
    attr(CI.PROB, 'conf.level')       <- conf.level
    attr(CI.PROB, 'method')           <- CI.method
    attr(CI.PROB, 'bootstrap.sample') <- MLE.prob.bootstrap
    attr(CI.PROB, 'bootstrap.sims')   <- bootstrap.sims
    
    #Add confidence interval to output
    OUT.MLE$CI.prob <- CI.PROB 
    
    } else {
      
      #Give warning if Hessian is small
      if (SE.prob/min(MLE.prob, 1-MLE.prob) > 10^2) {
        WARNING <- paste0('The computed Hessian of the log-likelihood is small ---\n', 
                          '    The assumptions for the asymptotic confidence interval may not hold \n', 
                          '    The asymptotic confidence interval is unreliable in this case \n', 
                          '    We recommend switching to bootstrap confidence interval by setting CI.method = \'bootstrap\'')
        warning(WARNING) } }
  
  #Set object class
  class(OUT.MLE) <- c('list', 'mle.matching')
  
  #Give output
  OUT.MLE }


print.mle.matching <- function(x, digits = 6, ...) {
  
  #Check input
  if (!('mle.matching' %in% class(x)))    stop('Error: This print method is for objects of class \'mle.matching\'')
  
  #Extract information
  DATA.NAME <- x$data.name
  DATA      <- x$data
  SIZE      <- x$size
  MLE.DF    <- x$MLE
  MLE.prob  <- MLE.DF[1, 1]
  MLE.phi   <- MLE.DF[1, 2]
  LL.MAX    <- x$maxloglike
  LL.MAX.M  <- x$maxloglike.mean
  CONF      <- x$CI.prob
  CI.method <- attributes(CONF)$method
  
  #Print title
  cat('\n    Maximum likelihood estimator (MLE) \n \n')
  cat(paste0('Data vector ', DATA.NAME, ' containing ', length(DATA), ' IID values from the generalised matching \n',
             'distribution with size = ', SIZE, ' and unknown probability parameter \n\n'))
  
  #Print MLE information
  cat(paste0('  MLE of probability parameter \n'))
  cat('    ', round(MLE.prob, digits), '\n\n')
  cat(paste0('  Maximum log-likelihood \n'))
  cat('    ', round(LL.MAX, digits), '\n\n')
  cat(paste0('  Likelihood-per-data-point (maximised geometric mean) \n'))
  cat('    ', round(exp(LL.MAX.M), digits), '\n\n')
  
  #Print confidence interval information
  if (CI.method == 'asymptotic') { cat('  Asymptotic') }
  if (CI.method == 'bootstrap')  { cat('  Bootstrap')  }
  cat(paste0(' ', 100*attributes(CONF)$conf.level, '% CI for probability parameter'))
  if (CI.method == 'bootstrap')  { 
    cat(paste0(' using ', attributes(CONF)$bootstrap.sims, ' resamples\n')) 
  } else {
    cat('\n') }
  cat(paste0('    [', round(CONF[1], digits), ', ', round(CONF[2], digits), ']\n\n')) }

