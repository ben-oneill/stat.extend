#' Optimal Confidence Intervals for finite populations
#' 
#' These functions compute an optimised confidence interval for statistics based on a sample.  The user may enter either a
#' data vector \code{x} or the sample size \code{n} and the sample statistic.  By default the confidence interval is computed
#' for an infinite population.  However, the user may enter a population size \code{N} and may use the logical value \code{unsampled} to specify
#' when the confidence interval is for the variance only of the unsampled part of the population.  This test accounts for the kurtosis, and
#' so the user must either specify the data vector or specify an assumed kurtosis \code{kurt}; if no kurtosis value is specified then the test
#' uses the sample kurtosis from the data.
#' 
#' The mean interval is built on a symmetric pivotal quantity so it is symmetric around the sample mean.  
#' 
#' The variance interval is built on a non-symmetric pivotal quantity, so it is optimised by taking the shortest possible confidence interval with the specified confidence level (see e.g., Tate and Klett 1959).  
#' 
#' The proportion interval uses the Wilson score interval (see e.g., Agresti and Coull 1998).
#' 
#' @param alpha alpha Numeric (probability) The significance level determining the confidence level for the interval (the confidence level is 1-alpha).
#' @param x  Numeric (vector) The vector of sample data. In the CONF.prop function this must be binary data. Ignored if a sample statistic is provided.
#' 
#' @param unsampled Logical (positive) Indicator of whether the user wants a confidence interval for the relevant parameter only for the unsampled part of the population (as opposed to the whole population)
#' @param kurt Numeric (positive) The assumed kurtosis of the underlying distribution (must be at least one)
#' 
#' @param n Integer (positive) The sample size
#' @param N Integer (positive) The population size (must be at least as large as the sample size)
#'
#' @return an object of class 'ci' providing the confidence interval and related information.
#'
#' @inheritParams checkIterArgs 
#' @examples 
#' DATA <- c(17.772, 16.359, 15.734, 15.698, 16.042, 
#' 15.527, 16.533, 15.385, 15.368, 18.603, 
#' 15.036, 13.873, 14.329, 15.837, 14.189, 
#' 15.398, 16.266, 12.970, 15.219, 16.444, 
#' 11.049, 14.262);
#' KURT <- 4.37559247659433 # moments::kurtosis(DATA);
#' CONF.mean(alpha = 0.1, x = DATA, N = 3200, kurt = KURT);
#' CONF.var(alpha = 0.1, x = DATA, N = 3200, kurt = KURT);
#' CONF.prop(alpha = 0.1, x = DATA > 15, N = 3200);
#' 
#' @name CONF
NULL


#'@rdname CONF
#' @param sample.mean Numeric (any) The sample mean of the data.
CONF.mean <- function(alpha, x = NULL, sample.mean = mean(x), 
                      sample.variance = var(x), n = length(x),
                      N = Inf, kurt = 3, unsampled = FALSE,
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check input alpha
  if (!is.numeric(alpha))   { stop('Error: alpha should be numeric') }
  if (length(alpha) != 1)   { stop('Error: alpha should be a single value'); }
  if (alpha < 0)            { stop('Error: alpha is negative'); }
  if (alpha > 1)            { stop('Error: alpha is greater than one'); }
  
  #Check congruence of data inputs
  if ((!missing(x) && !missing(sample.mean))) {
    if (abs(sample.mean - mean(x)) < 1e-15) 
      warning('specify data or sample mean but not both') else
        stop('Error: specify data or sample mean but not both'); }
  if ((!missing(x) && !missing(sample.variance))) {
    if (abs(sample.variance - var(x)) < 1e-15) 
      warning('specify data or sample variance but not both') else
        stop('Error: specify data or sample variance but not both'); }
  if ((!missing(x) && !missing(n))) {
    if (n != length(x)) 
      stop('Error: specify data or n but not both'); }
  
  #Check data inputs
  if (!is.numeric(x))     { stop('Error: x should be numeric') }
  if (!is.numeric(sample.mean))    { stop('Error: mean should be numeric') }
  if (length(sample.mean) != 1)    { stop('Error: mean should be a single value'); }
  if (!is.numeric(sample.variance)) { stop('Error: variance should be numeric') }
  if (length(sample.variance) != 1) { stop('Error: variance should be a single value'); }
  if (sample.variance < 0)          { stop('Error: variance is negative'); }
  if (!is.numeric(n))       { stop('Error: n should be numeric') }
  if (as.integer(n) != n)   { stop('Error: n should be an integer') }
  if (length(n) != 1)       { stop('Error: n should be a single value'); }
  if (n < 3)                { stop('Error: Method requires n > 2'); }
  
  #Check inputs N, kurt and unsampled
  if (!is.numeric(N))       { stop('Error: N should be numeric') }
  if (N != Inf) { 
    if (as.integer(N) != N) { stop('Error: N should be an integer') } }
  if (length(N) != 1)       { stop('Error: N should be a single value'); }
  if (N < n)                { stop('Error: N cannot be smaller than n'); }
  if (!is.numeric(kurt))    { stop('Error: kurt should be numeric') }
  if (length(kurt) != 1)    { stop('Error: kurt should be a single value'); }
  if (kurt < 1)             { stop('Error: kurt is less than one'); }
  if (!is.logical(unsampled)) { stop('Error: unsampled should be TRUE/FALSE') }
  if (length(unsampled) != 1) { stop('Error: unsampled should be a single value'); }
  
  #Check inputs for nlm
  checkIterArgs(gradtol, steptol, iterlim)

  #############
  
  #Compute CONF
  df  <- 2*n/(kurt - (n-3)/(n-1));
  TT  <- abs(qt(alpha/2, df, ncp = 0));
  UU  <- if (N == Inf) { 1 } else { (N-n)/N };
  ADJ <- if (!unsampled) { sqrt(UU) } else { 1/sqrt(UU) }
  L   <- sample.mean - sqrt(sample.variance)*ADJ*TT/sqrt(n);
  U   <- sample.mean + sqrt(sample.variance)*ADJ*TT/sqrt(n);
  CONF <- sets::interval(l = L, r = U, bounds = 'closed');
  attr(CONF, 'method') <- as.character(NA);
  
  #Add the description of the data
  if (is.null(x)) {
    DATADESC <- paste0('Interval uses ', n, ' data points with sample mean = ',
                       sprintf(sample.mean, fmt = '%#.4f'), 
                       ', sample variance ', 
                       sprintf(sample.variance, fmt = '%#.4f'),
                       ' and assumed kurtosis = ', 
                       sprintf(kurt, fmt = '%#.4f')) } else {
                         DATADESC <- paste0('Interval uses ', n, ' data points from data ', 
                                            deparse(substitute(x)), ' with sample variance = ', 
                                            sprintf(sample.variance, fmt = '%#.4f'),
                                            ' and assumed kurtosis = ', 
                                            sprintf(kurt, fmt = '%#.4f')); }
  attr(CONF, 'data') <- DATADESC;
  
  #Add class and attributes
  class(CONF) <- c('ci', 'interval');
  attr(CONF, 'confidence') <- 1 - alpha;
  attr(CONF, 'parameter')   <- if (N == Inf) {
    paste0('mean parameter for infinite population') } else {
      if (unsampled) {
        paste0('mean for unsampled population of size ', N-n) } else {
          paste0('mean for population of size ', N) } }
  
  CONF; }

#' @rdname CONF
#' @param sample.variance Numeric (non-neg) The sample variance of the data.
CONF.var <- function(alpha, x = NULL,
                     sample.variance = var(x), n = length(x),
                     N = Inf, kurt = NULL, unsampled = FALSE,
                     gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check input alpha
  if (!is.numeric(alpha))   { stop('Error: alpha should be numeric') }
  if (length(alpha) != 1)   { stop('Error: alpha should be a single value'); }
  if (alpha < 0)            { stop('Error: alpha is negative'); }
  if (alpha > 1)            { stop('Error: alpha is greater than one'); }
  
  #Check congruence of data inputs
  if ((!missing(x) && !missing(sample.variance))) {
    if (abs(sample.variance - var(x)) < 1e-15)
      warning('specify data or sample variance but not both') else
        stop('Error: specify data or sample variance but not both'); }
  if ((!missing(x) && !missing(n))) {
    if (n != length(x))
      stop('Error: specify data or n but not both'); }
  
  #Check data inputs
  if (!is.numeric(x))     { stop('Error: x should be numeric') }
  if (!is.numeric(sample.variance)) { stop('Error: variance should be numeric') }
  if (length(sample.variance) != 1) { stop('Error: variance should be a single value'); }
  if (sample.variance < 0)          { stop('Error: variance is negative'); }
  if (!is.numeric(n))       { stop('Error: n should be numeric') }
  if (as.integer(n) != n)   { stop('Error: n should be an integer') }
  if (length(n) != 1)       { stop('Error: n should be a single value'); }
  if (n < 3)                { stop('Error: Method requires n > 2'); }
  
  #Check inputs N, kurt and unsampled
  if (!is.numeric(N))       { stop('Error: N should be numeric') }
  if (N != Inf) {
    if (as.integer(N) != N) { stop('Error: N should be an integer') } }
  if (length(N) != 1)       { stop('Error: N should be a single value'); }
  if (N <= n)               { stop('Error: N should be larger than n'); }
  miss.kurt <- missing(kurt);
  if(!miss.kurt)              {
    if (!is.numeric(kurt))    { stop('Error: kurt should be numeric') }
    if (length(kurt) != 1)    { stop('Error: kurt should be a single value'); }
    if (kurt < 1)             { stop('Error: kurt is less than one'); } }
  if (!is.logical(unsampled)) { stop('Error: unsampled should be TRUE/FALSE') }
  if (length(unsampled) != 1) { stop('Error: unsampled should be a single value'); }
  if ((unsampled) & (N-n < 3)) {
    stop('Error: Method requires N-n > 2'); }
  
  #Check inputs for nlm
  if (!is.numeric(gradtol)) { stop('Error: gradtol should be numeric') }
  if (length(gradtol) != 1) { stop('Error: gradtol should be a single value'); }
  if (gradtol <= 0)         { stop('Error: gradtol should be positive'); }
  if (!is.numeric(steptol)) { stop('Error: steptol should be numeric') }
  if (length(steptol) != 1) { stop('Error: steptol should be a single value'); }
  if (steptol <= 0)         { stop('Error: steptol should be positive'); }
  if (!is.numeric(iterlim)) { stop('Error: iterlim should be numeric') }
  if (length(iterlim) != 1) { stop('Error: iterlim should be a single value'); }
  if (iterlim <= 0)         { stop('Error: iterlim should be positive'); }
  
  #############
  
  #Determine the kurtosis value
  if (miss.kurt) {
    if (missing(x)) { kurt <- 3; } else {
      kurt <- n*sum((x-mean(x))^4)/(sum((x-mean(x))^2)^2) } }
  
  #Compute the confidence interval (including sample data)
  if (!unsampled) {
    
    #Simplify probability functions (with stipulated parameters)
    df1 <- 2*n/(kurt - (n-3)/(n-1));
    df2 <- 2*(N-n)/(2 + (kurt-3)*(1-2/N-1/(N*n)));
    Q <- function(L) { qf(L, df1, df2); }
    f <- function(L) { df(L, df1, df2); }
    
    #Set objective function
    WW <- function(phi) {
      
      #Set parameter functions
      T0 <- alpha/(1+exp(-phi));
      T1 <- T0/(1+exp(phi));
      T2 <- T1*((1-exp(2*phi))/(1+2*exp(phi)+exp(2*phi)));
      
      #Set interval bounds and objective
      L  <- Q(T0);
      U  <- Q(T0 + 1 - alpha);
      W0 <- 1/L - 1/U;
      
      #Set gradient of objective
      if (!is.null(f)) {
        attr(W0, 'gradient') <- T1*(1/(f(U)*U^2) - 1/(f(L)*L^2)); }
      
      W0; }
    
    #Compute the HDR
    #The starting value for the parameter phi is set to zero
    #This is the exact optima in the case of a symmetric distribution
    OPT <- nlm(WW, p = 0,
               gradtol = gradtol, steptol = steptol, iterlim = iterlim);
    TT <- alpha/(1+exp(-OPT$estimate));
    A <- (n-1)/(N-1);
    B <- if (N == Inf) { 1 } else { (N-n)/(N-1) }
    L   <- A + B/Q(TT + 1 - alpha);
    U   <- A + B/Q(TT);
    CONF <- sample.variance*sets::interval(l = L, r = U, bounds = 'closed');
    
    #Add the description of the method
    METHOD <- ifelse((OPT$iterations == 1),
                     paste0('Computed using nlm optimisation with ',
                            OPT$iterations, ' iteration (code = ', OPT$code, ')'),
                     paste0('Computed using nlm optimisation with ',
                            OPT$iterations, ' iterations (code = ', OPT$code, ')'));
    attr(CONF, 'method') <- METHOD; }
  
  #Compute the confidence interval (excluding sample data)
  if (unsampled) {
    
    #Simplify probability functions (with stipulated parameters)
    df1 <- 2*n/(kurt - (n-3)/(n-1));
    df2 <- 2*(N-n)/(kurt - (N-n-3)/(N-n-1));
    Q <- function(L) { qf(L, df1, df2); }
    f <- function(L) { df(L, df1, df2); }
    
    #Set objective function
    WW <- function(phi) {
      
      #Set parameter functions
      T0 <- alpha/(1+exp(-phi));
      T1 <- T0/(1+exp(phi));
      T2 <- T1*((1-exp(2*phi))/(1+2*exp(phi)+exp(2*phi)));
      
      #Set interval bounds and objective
      L  <- Q(T0);
      U  <- Q(T0 + 1 - alpha);
      W0 <- 1/L - 1/U;
      
      #Set gradient of objective
      if (!is.null(f)) {
        attr(W0, 'gradient') <- T1*(1/(f(U)*U^2) - 1/(f(L)*L^2)); }
      
      W0; }
    
    #Compute the HDR
    #The starting value for the parameter phi is set to zero
    #This is the exact optima in the case of a symmetric distribution
    OPT <- nlm(WW, p = 0,
               gradtol = gradtol, steptol = steptol, iterlim = iterlim);
    TT <- alpha/(1+exp(-OPT$estimate));
    L   <- 1/Q(TT + 1 - alpha);
    U   <- 1/Q(TT);
    CONF <- sample.variance*sets::interval(l = L, r = U, bounds = 'closed');
    
    #Add the description of the method
    METHOD <- ifelse((OPT$iterations == 1),
                     paste0('Computed using nlm optimisation with ',
                            OPT$iterations, ' iteration (code = ', OPT$code, ')'),
                     paste0('Computed using nlm optimisation with ',
                            OPT$iterations, ' iterations (code = ', OPT$code, ')'));
    attr(CONF, 'method') <- METHOD; }
  
  #Add the description of the data
  if (is.null(x))  {
    DATADESC <- paste0('Interval uses ', n,
                       ' data points with sample variance = ',
                       sprintf(sample.variance, fmt = '%#.4f'),
                       ' and assumed kurtosis = ',
                       sprintf(kurt, fmt = '%#.4f')) }
  if (!is.null(x)) {
    if (miss.kurt) {
      DATADESC <- paste0('Interval uses ', n, ' data points from data ',
                         deparse(substitute(x)), ' with sample variance = ',
                         sprintf(sample.variance, fmt = '%#.4f'),
                         ' and sample kurtosis = ',
                         sprintf(kurt, fmt = '%#.4f')); } else {
                           DATADESC <- paste0('Interval uses ', n, ' data points from data ',
                                              deparse(substitute(x)), ' with sample variance = ',
                                              sprintf(sample.variance, fmt = '%#.4f'),
                                              ' and assumed kurtosis = ',
                                              sprintf(kurt, fmt = '%#.4f')); } }
  attr(CONF, 'data') <- DATADESC;
  
  #Add class and attributes
  class(CONF) <- c('ci', 'interval');
  attr(CONF, 'confidence')  <- 1 - alpha;
  attr(CONF, 'parameter')   <- if (N == Inf) {
    paste0('variance parameter for infinite population') } else {
      if (unsampled) {
        paste0('variance for unsampled population of size ', N-n) } else {
          paste0('variance for population of size ', N) } }
  
  CONF; }






#' @param sample.prop Numeric (probability) The sample proportion of the data (only for binary data)
#' @rdname CONF
CONF.prop <- function(alpha, x = NULL, 
                      sample.prop = mean(x), n = length(x),
                      N = Inf, unsampled = FALSE) {
  
  #Check input alpha
  if (!is.numeric(alpha))   { stop('Error: alpha should be numeric') }
  if (length(alpha) != 1)   { stop('Error: alpha should be a single value'); }
  if (alpha < 0)            { stop('Error: alpha is negative'); }
  if (alpha > 1)            { stop('Error: alpha is greater than one'); }
  
  #Check congruence of data inputs
  if ((!missing(x) && !missing(sample.prop))) {
    if (abs(sample.prop - mean(x)) < 1e-15) 
      warning('specify data or sample proportion but not both') else
        stop('Error: specify data or sample proportion but not both'); }
  if ((!missing(x) && !missing(n))) {
    if (n != length(x)) 
      stop('Error: specify data or n but not both'); }
  
  #Check data inputs
  P <- sample.prop;
  if (!missing(x)) {
    xexpr <- deparse(substitute(x));
    if (is.logical(x)) x <- x + 0;
    if (!is.numeric(x))     { stop('Error: x should be numeric') }
    if (!all(x %in% c(0L,1L))) { stop('Error: x should be binary data') } }
  if (!is.numeric(P))       { stop('Error: sample proportion should be numeric') }
  if (length(P) != 1)       { stop('Error: sample proportion should be a single value'); }
  if (!is.numeric(n))       { stop('Error: n should be numeric') }
  if (as.integer(n) != n)   { stop('Error: n should be an integer') }
  if (length(n) != 1)       { stop('Error: n should be a single value'); }
  if (n < 1)                { stop('Error: Must have at least one data point'); }
  
  #Compute the confidence interval in trivial cases
  if (alpha == 0) {
    CONF <- sets::interval(l = 0, r = 1, bounds = 'closed'); }
  if (alpha == 1) {
    CONF <- sets::interval(l = P, r = P, bounds = 'closed'); }
  
  #Compute the confidence interval in non-trivial case
  if ((alpha > 0) && (alpha < 1)) {
    if (N == Inf) {
      psi <- qchisq(1-alpha, df = 1)/(2*n) } else {
        if (unsampled) {
          psi <- ((N-1)/(N-n))*qchisq(1-alpha, df = 1)/(2*n) } else {
            psi <- ((N-n)/(N-1))*qchisq(1-alpha, df = 1)/(2*n) } }
    T1  <- (P+psi)/(1+2*psi);
    T2  <- sqrt(2*psi*P*(1-P) + psi^2)/(1+2*psi);
    L   <- T1 - T2;
    U   <- T1 + T2;
    CONF <- sets::interval(l = L, r = U, bounds = 'closed'); }
  
  #Add the description of the data
  if (missing(x)) {
    DATADESC <- paste0('Interval uses ', n, 
                       ' binary data points with sample proportion = ', 
                       sprintf(P, fmt = '%#.4f')) } else {
                         DATADESC <- paste0('Interval uses ', n, ' binary data points from data ', 
                                            xexpr, ' with sample proportion = ', 
                                            sprintf(P, fmt = '%#.4f')) }
  attr(CONF, 'data') <- DATADESC;
  
  #Add class and attributes
  class(CONF) <- c('ci', 'interval');
  attr(CONF, 'method')      <- as.character(NA);
  attr(CONF, 'confidence')  <- 1 - alpha;
  attr(CONF, 'parameter')   <- if (N == Inf) {
    paste0('proportion parameter for infinite population') } else {
      if (unsampled) {
        paste0('proportion for unsampled population of size ', N-n) } else {
          paste0('proportion for population of size ', N) } }
  
  CONF; }

print.ci <- function(x, ...) {
  
  #Print description of confidence interval
  cat('\n        Confidence Interval (CI) \n \n');
  cat(paste0(sprintf(100*attributes(x)$confidence, fmt = '%#.2f'), '%'),
      'CI for', attributes(x)$parameter, '\n');
  
  #Print data description
  cat(attributes(x)$data, '\n');
  
  #Print method
  if (!is.na(attributes(x)$method)) {
    cat(attributes(x)$method, '\n'); }
  
  #Print confidence interval
  cat('\n');
  writeLines(as.character(c(x)))
  invisible(c(x))
  cat('\n'); }
