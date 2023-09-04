######
## Penalized.cor.gof -- The Penalized ACF/PACF estimator for time series goodness of fit.
##
## Introduction for the arguments:
##
##     x        -- a numeric vector or univariate time series
##     lag      -- the lag which the statistic will be based on to calculate penalized acf/pacf, default = 10
##     type     -- the type of test which penalized acf/pacf to be used, default = "LB"
##    fitdf     -- the number of degrees of freedom to be subtracted if x is a series of residuals, default = 0
##    alpha     -- the nominal level selected, default = 0.05
##
## The following are types of acf/pacf used in tests
##     BP -- Box Pierce test; WBP -- Weighted Box Pierce test;
##     LB -- Ljung Box test; WLB -- Weighted Ljung Box test; M -- Monti test; WM -- Weighted Monti test;
##     MM -- Mahdi McLeod test; LM -- Li-Mak test; WLM -- Weighted Li-Mak test;
##     LM.pacf -- Li-Mak test with pacf; WLM.pacf -- Weighted Li-Mak test with pacf
## Note that we have lag-fitdf penalized acf/pacfs if type = "LM","WLM","LM.pacf","WLM.pacf".
######

Penalized.cor.gof <- function(x, lag = 10, fitdf= 0,
                              type = c("LB","BP","M","WLB","WBP","WM","MM","LM","WLM","LM.pacf","WLM.pacf"),
                              alpha = 0.05){
  if(!is.numeric(x))
    stop("'x' must be numeric")
  if (NCOL(x) > 1)
    stop("'x' is not a vector or univariate time series")
  if (lag < 1)
    stop("'lag' must be positive")
  if (fitdf < 0)
    stop("'fitdf' cannot be negative")
  if (fitdf >= lag)
    stop("'lag' must exceed fitted degrees of freedom 'fitdf'")
  
  series <- deparse(substitute(x))
  n <- as.integer(sum(!is.na(x)))
  m <- as.integer(lag)
  
  ## find the type and obtain the acf/pacf
  type <- match.arg(type)
  if (type == "LM" | type =="WLM" | type=="LM.pacf" | type =="WLM.pacf"){
    k <- (fitdf + 1) : m
  }else{
    k <- 1 : m
  }
  if(type == "M" | type == "WM" | type == "MM" | type =="LM.pacf" | type =="WLM.pacf"){
    rhat <- stats::pacf(x, lag.max=m, plot=FALSE, na.action=na.pass)$acf[1:m]
    rhat <- rhat[k]
    r <- abs(rhat)
  }else{
    rhat <- stats::acf(x, lag.max=m, plot=FALSE, na.action=na.pass)$acf[2:(m+1)]
    rhat <- rhat[k]
    r <- abs(rhat)
  }
  
  if (type == "BP"){
    test <- "Box-Pierce"
    wk <- 1
    qalpha <- qchisq(1-alpha,df = m-fitdf)
    lk <- sqrt(1/n)
    uk <- sqrt(qalpha/(n*wk))
  }
  if (type == "WBP"){
    test <- "Weighted Box-Pierce"
    wk <- (m + 1 - k)/m
    shape <- (3/4) * (m + 1)^2 * m/(2 * m^2 + 3 * m + 1 - 6 * m * fitdf)
    scale <- (2/3) * (2 * m^2 + 3 * m + 1 - 6 * m * fitdf)/m/(m + 1)
    qalpha <- qgamma(1-alpha, shape = shape, scale = scale)
    lk <- sqrt(1/n)
    u <- (qalpha/m)
    uk <- sqrt(u/(n*wk))
  }
  if (type == "LB"){
    test <- "Ljung-Box"
    wk <- (n + 2)/((n - k))
    qalpha <- qchisq(1-alpha,df = m-fitdf)
    lk <- sqrt((n-k)/(n^2+2*n))
    u <- (qalpha/m)
    uk <- sqrt(u/(n*wk))
  }
  if(type == "WLB"){
    test <- "Weighted Ljung-Box"
    wk <- (n + 2)*(m+1-k)/(m*(n - k))
    lk <- sqrt((n - k)/(n^2 + 2*n))
    shape <- (3/4) * (m + 1)^2 * m/(2 * m^2 + 3 * m + 1 - 6 * m * fitdf)
    scale <- (2/3) * (2 * m^2 + 3 * m + 1 - 6 * m * fitdf)/m/(m + 1)
    qalpha <- qgamma(1-alpha, shape = shape, scale = scale)
    u <- (qalpha/m)
    uk <- sqrt(u/(n*wk))
    lk <- pmin(lk,uk-.001)
  }
  if(type == "M"){
    test <- "Monti"
    wk <-  (n + 2)/(n - k)
    qalpha <- qchisq(1-alpha, df = m-fitdf)
    u <- qalpha/m
    lk <- sqrt((n - k)/(n^2 + 2*n))
    uk <- sqrt(u/(n*wk))
  }
  if (type == "WM"){
    test <- "Weighted Monti"
    wk <- (n + 2)*(m + 1 - k)/(m*(n - k))
    lk <- sqrt((n - k)/(n^2 + 2*n))
    shape <- (3/4) * (m + 1)^2 * m/(2 * m^2 + 3 * m + 1 - 6 * m * fitdf)
    scale <- (2/3) * (2 * m^2 + 3 * m + 1 - 6 * m * fitdf)/m/(m + 1)
    qalpha <- qgamma(1-alpha, shape = shape, scale = scale)
    u <- qalpha/m
    uk <- sqrt(u/(n*wk))
    lk <- pmin(lk,uk-.001)
  }
  if(type == "MM"){
    test <- "Mahdi-McLeod"
    wk <- 3*(m + 1 - k) /(2*m + 1)  
    df.mm <- 1.5 * m * (m + 1)/(2*m + 1) - fitdf
    qalpha <- qchisq(1-alpha,df = df.mm)
    u <- (qalpha/m)
    lk <- sqrt((n - k)/(n^2 + 2*n))
    uk <- sqrt(1 - exp(-qalpha*(2*m + 1)/(3*n*m*(m + 1 - k))))
    
    lk <- pmin(lk,uk-.001)
  }
  if(type == "LM"){
    test <- "Li-Mak on autocorrelations"
    wk <- 1
    qalpha <- qchisq(1-alpha,df = m-fitdf)
    lk <- sqrt(1/n)
    uk <- sqrt(qalpha/(m - fitdf))*lk
  }
  if(type == "WLM"){
    test <- "Weighted Li-Mak on autocorrelations"
    wk <- (m - k + fitdf + 1)/m
    shape <- (3/4) * (m + fitdf + 1)^2 * (m - fitdf)/(2 * m^2 + 3 * m +
                                                        2 * m * fitdf + 2 * fitdf^2 +
                                                        3 * fitdf + 1)
    scale <- (2/3) * (2 * m^2 + 3 * m + 2 * m * fitdf +
                        2 * fitdf^2 + 3 * fitdf + 1)/(m * (m + fitdf + 1))
    qalpha <- qchisq(1-alpha,df = m-fitdf)
    lk <- sqrt(1/n)
    l2  <- sqrt(1/(n*(m - k + (fitdf + 1))/m))
    uk <- sqrt(qalpha/(m-fitdf))*l2
  }
  if(type == "LM.pacf"){
    test <- "Li-Mak on partial autocorrelations"
    wk <- 1
    qalpha <- qchisq(1-alpha,df = m-fitdf)
    lk <- sqrt(1/n)
    uk <- sqrt(qalpha/(m - fitdf))*lk
  }
  if(type == "WLM.pacf"){
    test <- "Weighted Li-Mak on partial autocorrelations"
    wk <- (m - k + fitdf + 1)/m
    shape <- (3/4) * (m + fitdf + 1)^2 * (m - fitdf)/(2 * m^2 + 3 * m +
                                                        2 * m * fitdf + 2 * fitdf^2 +
                                                        3 * fitdf + 1)
    scale <- (2/3) * (2 * m^2 + 3 * m + 2 * m * fitdf +
                        2 * fitdf^2 + 3 * fitdf + 1)/(m * (m + fitdf + 1))
    qalpha <- qgamma(1-alpha, shape = shape, scale = scale)
    qalpha <- qchisq(1-alpha,df = m-fitdf)
    lk <- sqrt(1/n)
    uk <- sqrt(qalpha/(n*wk*(m-fitdf)))
  }
  rtar <- (1 + (r - lk)/lk) * r  * (sign(lk - r) + 1)/2 + 
    (r + (r - lk)*(uk - r)/((uk - lk)))*(sign(r - lk) + 1)/2*(sign(uk - r) + 1)/2 +
    r*(1 + (m-fitdf)/(m*(1-alpha))*sqrt((m + 1 - k)/(wk*m))*sqrt(n)/(2*log(n))*(r - uk)*(1 - r)/(1 - uk))*(sign(r - uk) + 1)/2
  rtar <- sign(rhat) * rtar
  lambda <- 0
  if (m == 1){
    lambda <- 0
  }else{
    lambda <- m/(4*log(m)*wk)*(1 + abs(r-lk)/r^2)* (sign(lk - r) + 1)/2 +
      (1 + sqrt(n)/(4*log(m))*abs(r - uk)*abs(r-lk)/abs(uk - lk)^2)*(sign(r - lk) + 1)/2*(sign(uk - r) + 1)/2 + 
      sqrt(n)/(4*log(m))*(m + 1 - k)/(m*wk)*(1 + (r - lk)^2/((r - uk)*(1 - uk)))*(sign(r - uk) + 1)/2
  }
  w <- lambda/(1 + lambda)
  rtilde <- w*rtar + (1 - w)*rhat
  structure(list(r.pen = rtilde, type = test, n.used = n,
                 lag = lag, fitdf = fitdf, series.name = series))
}


######
## Penalized.Box.test -- The Penalized (Weighted) Box-Pierce, (Weighted) Ljung-Box, (Weighted) Monti, Weighted McLeod-Li type test for fitted ARMA and detection of nonlinear Models
##
## Introduction for the arguments:
##
##      x       -- the residuals or initial data which is a numeric vector or univariate time series
##     lag      -- the lag which the statistic will be based on to calculate penalized acf/pacf, default = 10
##     type     -- the type of test, default = "Ljung-Box"
##    fitdf     -- the number of degrees of freedom to be subtracted if x is a series of residuals, default = 0
##    alpha     -- the nominal level selected, default = 0.05
##
## The followings are used to perform transformations to detect nonlinear processes, all default=FALSE
##
##    squared  -- Take the squared value of residuals.
##    log.squared -- Take the log of the squared residuals
##    absolute   -- Take the absolute value of the residuals
##
##  For backward compatability
##
##   penalized  -- If FALSE, perform the unpenalized tests using the sample ACF/PACF estimates, default=TRUE
##   weighted   -- If FALSE, perform the original Box-Pierce, Ljung-Box or Monti or McLeod Li type if using a transformation, default=TRUE.
######

Penalized.Box.test <- function(x, lag = 10, type = c("Ljung-Box", "Box-Pierce", "Monti"),
                                        fitdf = 0, squared = FALSE, log.squared = FALSE, absolute = FALSE,
                                        weighted = TRUE, penalized = TRUE, alpha=0.05) {
  if(!is.numeric(x))
    stop("'x' must be numeric")
  if (NCOL(x) > 1)
    stop("'x' is not a vector or univariate time series")
  if (lag < 1)
    stop("'lag' must be positive")
  if (fitdf < 0)
    stop("'fitdf' cannot be negative")
  if (fitdf >= lag)
    stop("'lag' must exceed fitted degrees of freedom 'fitdf'")
  if ((squared && log.squared) || (squared && absolute) ||
      (log.squared && absolute))
    stop("only one option of 'squared', 'log.squared' or 'absolute' can be selected")

  DNAME <- deparse(substitute(x))
  type <- match.arg(type)
  n <- as.integer(sum(!is.na(x)))
  lag <- as.integer(lag)

  ## transform residuals for detection of nonlinear models
  if (squared) {
    x <- x^2
    fitdf <- 0
  }
  if (log.squared) {
    x <- log(x^2)
    fitdf <- 0
  }
  if (absolute) {
    x <- abs(x)
    fitdf <- 0
  }

  if (weighted) { #Weighted
    if (type == "Box-Pierce"){
      if (penalized){  #Penalized
        METHOD <- "Penalized weighted Box-Pierce test (Gamma Approximation)"
        cor <- Penalized.cor.gof(x, lag=lag, fitdf=fitdf, type="WBP", alpha=alpha)
        obs <- cor$r.pen[1:lag]
      }
      else{  #Not penalized
        METHOD <- "Weighted Box-Pierce test (Gamma Approximation)"
        cor <- stats::acf(x, lag.max = lag, plot=FALSE, na.action=na.pass)
        obs <- cor$acf[2:(lag + 1)];
      }
    }
    if (type =="Ljung-Box"){
      if (penalized){
        METHOD <- "Penalized weighted Ljung-Box test (Gamma Approximation)"
                cor <- Penalized.cor.gof(x, lag=lag, fitdf=fitdf, type="WLB", alpha=alpha)
        obs <- cor$r.pen[1:lag]
      }
      else{
        METHOD <- "Weighted Ljung-Box test (Gamma Approximation)"
        cor <- stats::acf(x, lag.max = lag, plot=FALSE, na.action=na.pass)
        obs <- cor$acf[2:(lag + 1)]
      }
    }
    if (type == "Monti") {
      if (penalized){
        METHOD <- "Penalized weighted Monti test (Gamma Approximation)"
        cor <- Penalized.cor.gof(x, lag=lag, fitdf=fitdf, type="WM", alpha=alpha)
        obs <- cor$r.pen[1:lag]
      }
      else{
        METHOD <- "Weighted Monti test (Gamma Approximation)"
        cor <- stats::pacf(x, lag.max = lag,plot=FALSE, na.action=na.pass)
        obs <- cor$acf[1:lag]
      }
    }

    weights <- (lag - 1:lag + 1)/(lag)
    if (type == "Box-Pierce") {
      STATISTIC <- n * sum(weights * obs^2)
    }
    else {
      STATISTIC <- n * (n + 2) * sum(weights * (1/seq.int(n - 1, n - lag) * obs^2))
    }

    if (penalized){
      if (squared) {
        names(STATISTIC) <- "Penalized weighted X-squared on Squared Residuals for detecting nonlinear processes"
      }
      else if (log.squared) {
        names(STATISTIC) <- "Penalized weighted X-squared on Log-squared Residuals for detecting nonlinear processes"
      }
      else if (absolute) {
        names(STATISTIC) <- "Penalized weighted X-squared on Absolute valued Residuals for detecting nonlinear processes"
      }
      else {
        names(STATISTIC) <- "Penalized weighted X-squared on Residuals for fitted ARMA process"
      }
    }else{
      if (squared) {
        names(STATISTIC) <- "Weighted X-squared on Squared Residuals for detecting nonlinear processes"
      }
      else if (log.squared) {
        names(STATISTIC) <- "Weighted X-squared on Log-squared Residuals for detecting nonlinear processes"
      }
      else if (absolute) {
        names(STATISTIC) <- "Weighted X-squared on Absolute valued Residuals for detecting nonlinear processes"
      }
      else {
        names(STATISTIC) <- "Weighted X-squared on Residuals for fitted ARMA process"
      }
    }

    shape <- (3/4) * (lag + 1)^2 * lag/(2 * lag^2 + 3 * lag + 1 - 6 * lag * fitdf)
    scale <- (2/3) * (2 * lag^2 + 3 * lag + 1 - 6 * lag * fitdf)/lag/(lag + 1)
    PARAMETER <- c(Shape=shape, Scale=scale)
    PVAL <- 1 - pgamma(STATISTIC, shape = shape, scale = scale)
    names(PVAL) <- "Approximate p-value"
  }
  else { #Not weighted
    if (type == "Box-Pierce"){
      if (penalized){ #Penalized
        METHOD <- "Penalized Box-Pierce test"
        cor <- Penalized.cor.gof(x, lag=lag, fitdf=fitdf, type="BP", alpha=alpha)
        obs <- cor$r.pen[1:lag]
      }
      else{  #Not Penalized
        METHOD <- "Box-Pierce test"
        cor <- stats::acf(x, lag.max = lag, plot=FALSE, na.action=na.pass)
        obs <- cor$acf[2:(lag + 1)]
      }
    }
    if (type =="Ljung-Box"){
      if (penalized){
        METHOD <- "Penalized Ljung-Box test"
        cor <- Penalized.cor.gof(x, lag=lag, fitdf=fitdf, type="LB", alpha=alpha)
        obs <- cor$r.pen[1:lag]
      }
      else{
        METHOD <- "Ljung-Box test"
        cor <- stats::acf(x, lag.max = lag, plot=FALSE, na.action=na.pass)
        obs <- cor$acf[2:(lag + 1)]
      }
    }
    if (type == "Monti") {
      if (penalized){
        METHOD <- "Penalized Monti test"
        cor <- Penalized.cor.gof(x, lag=lag, fitdf=fitdf, type="M", alpha=alpha)
        obs <- cor$r.pen[1:lag]
      }
      else{
        METHOD <- "Monti test"
        cor <- stats::pacf(x, lag.max = lag, plot=FALSE, na.action=na.pass)
        obs <- cor$acf[1:lag];
      }
    }

    if (type == "Box-Pierce") {
      STATISTIC <- n * sum(obs^2)
    }
    else {
      STATISTIC <- n * (n + 2) * sum((1/seq.int(n - 1, n - lag) * obs^2))
    }

    if (penalized){
      if (squared) {
        names(STATISTIC) <- "Penalized X-squared on Squared Residuals for detecting nonlinear processes"
      }
      else if (log.squared) {
        names(STATISTIC) <- "Penalized X-squared on Log-squared Residuals for detecting nonlinear processes"
      }
      else if (absolute) {
        names(STATISTIC) <- "Penalized X-squared on Absolute valued Residuals for detecting nonlinear processes"
      }
      else {
        names(STATISTIC) <- "Penalized X-squared on Residuals for fitted ARMA process"
      }
    }else{
      if (squared) {
        names(STATISTIC) <- "X-squared on Squared Residuals for detecting nonlinear processes"
      }
      else if (log.squared) {
        names(STATISTIC) <- "X-squared on Log-squared Residuals for detecting nonlinear processes"
      }
      else if (absolute) {
        names(STATISTIC) <- "X-squared on Absolute valued Residuals for detecting nonlinear processes"
      }
      else {
        names(STATISTIC) <- "X-squared on Residuals for fitted ARMA process"
      }
    }


    PARAMETER <- c(df = lag - fitdf)
    PVAL <- 1 - pchisq(STATISTIC, lag - fitdf)
    names(PVAL) <- "p-value"
  }

  structure(list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = PVAL,
                 method = METHOD,
                 data.name = DNAME),
            class = "htest")
}

######
## Penalized.MahdiMcLeod.test -- The Penalized MahdiMcLeod test for fitted ARMA and detection of nonlinear Models
##
## Introduction for the arguments:
##
##      x       -- the residuals or initial data which is a numeric vector or univariate time series
##     lag      -- the lag which the statistic will be based on to calculate penalized acf/pacf, default = 10
##    fitdf     -- the number of degrees of freedom to be subtracted if x is a series of residuals, default = 0
##    alpha     -- the nominal level selected, default = 0.05
##
## The followings are used to perform transformations to detect nonlinear processes, all default=FALSE
##
##    squared  -- Take the squared value of residuals.
##    log.squared -- Take the log of the squared residuals
##    absolute   -- Take the absolute value of the residuals
##
## For backward compatability
##   penalized  -- If FALSE, perform the orginal Mahdi-McLeod test using the sample PACF estimates
######

Penalized.MahdiMcLeod.test <- function(x, lag = 10,  fitdf = 0, squared = FALSE, log.squared = FALSE,
                                       absolute = FALSE, penalized = TRUE, alpha=0.05){
  if(!is.numeric(x))
    stop("'x' must be numeric")
  if (NCOL(x) > 1)
    stop("'x' is not a vector or univariate time series")
  if (lag < 1)
    stop("'lag' must be positive")
  if (fitdf < 0)
    stop("'fitdf' cannot be negative")
  if (fitdf >= lag)
    stop("'lag' must exceed fitted degrees of freedom 'fitdf'")
  if ((squared && log.squared) || (squared && absolute) ||
      (log.squared && absolute))
    stop("only one option of 'squared', 'log.squared' or 'absolute' can be selected")


  DNAME <- deparse(substitute(x))
  n <- as.integer(sum(!is.na(x)))
  lag <- as.integer(lag)

  ## transform residuals for detection of nonlinear models
  if (squared) {
    x <- x^2
    fitdf <- 0
  }
  if (log.squared) {
    x <- log(x^2)
    fitdf <- 0
  }
  if (absolute) {
    x <- abs(x)
    fitdf <- 0
  }

  k <- 1 : lag
  if (penalized){
    METHOD <- "Penalized Mahdi-McLeod test"
    cor <- Penalized.cor.gof(x, lag=lag, fitdf=fitdf, type="MM", alpha=alpha)
    obs <- cor$r.pen[1:lag]
  }
  else{
    METHOD <- "Mahdi-McLeod test"
    cor <- stats::pacf(x, lag.max = lag, plot=FALSE, na.action=na.pass)
    obs <- cor$acf[1:lag]
  }

  STATISTIC <- (-3 * n/(2 * lag + 1)) * sum((lag + 1 - k)*log(1 - obs^2))

  if (penalized){
    if (squared) {
      names(STATISTIC) <- "Penalized X-squared on Squared Residuals for detecting nonlinear processes"
    }
    else if (log.squared) {
      names(STATISTIC) <- "Penalized X-squared on Log-squared Residuals for detecting nonlinear processes"
    }
    else if (absolute) {
      names(STATISTIC) <- "Penalized X-squared on Absolute valued Residuals for detecting nonlinear processes"
    }
    else {
      names(STATISTIC) <- "Penalized X-squared on Residuals for fitted ARMA process"
    }
  }else{
    if (squared) {
      names(STATISTIC) <- "X-squared on Squared Residuals for detecting nonlinear processes"
    }
    else if (log.squared) {
      names(STATISTIC) <- "X-squared on Log-squared Residuals for detecting nonlinear processes"
    }
    else if (absolute) {
      names(STATISTIC) <- "X-squared on Absolute valued Residuals for detecting nonlinear processes"
    }
    else {
      names(STATISTIC) <- "X-squared on Residuals for fitted ARMA process"
    }
  }

  mydf <- (1.5 * lag * (lag + 1)/(2 * lag + 1) - fitdf)
  PARAMETER <- c(df = mydf)
  PVAL <- 1 - pchisq(STATISTIC, mydf)
  names(PVAL) <- "p-value"
  structure(list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = PVAL,
                 method = METHOD,
                 data.name = DNAME),
            class = "htest")
}


######
## Penalized.LM.test -- The Penalized Li-Mak type test
##
## Introduction for the arguments:
##
##      x       -- the residuals or initial data which is a numeric vector or univariate time series
##     h.t      -- the sample conditional variances
##     lag      -- the lag which the statistic will be based on to calculate penalized acf/pacf, default = 10
##     type     -- The type of test based on acf or pacf, default="correlation"
##    fitdf     -- the number of degrees of freedom to be subtracted if x is a series of residuals, default = 0
##    alpha     -- the nominal level selected, default = 0.05
##
## For backward compatability
##
##   penalized  -- If FALSE, perform the unpenalized tests using the sample ACF/PACF estimates, default=TRUE
##   weighted   -- If FALSE, perform the Li-Mak test, default=TRUE.
######

Penalized.LM.test <- function (x, h.t, lag = 10, type = c("correlation", "partial"),
                                        fitdf = 0, weighted = TRUE, penalized = TRUE, alpha=0.05){
  if(!is.numeric(x))
    stop("'x' must be numeric")
  if (NCOL(x) > 1)
    stop("'x' is not a vector or univariate time series")
  if (lag < 1)
    stop("'lag' must be positive")
  if (fitdf < 0)
    stop("'fitdf' cannot be negative")
  if (fitdf >= lag)
    stop("'lag' must exceed fitted degrees of freedom 'fitdf'")

  DNAME <- deparse(substitute(x))
  type <- match.arg(type)
  n <- as.integer(sum(!is.na(x)))
  lag <- as.integer(lag)
  x <- x^2/h.t
  k <- (fitdf + 1):lag

  if (weighted) {
    if (type == "correlation"){
      if(penalized){
        METHOD <- "Penalized weighted Li-Mak test on autocorrelations (Gamma Approximation)"
        cor <- Penalized.cor.gof(x, lag=lag, fitdf=fitdf, type="WLM", alpha=alpha)
        obs <- cor$r.pen
      }
      else{
        METHOD <- "Weighted Li-Mak test on autocorrelations (Gamma Approximation)"
        cor <- stats::acf(x, lag.max = lag, type="correlation", plot=FALSE, na.action=na.pass)
        obs <- cor$acf[2:(lag + 1)][k]
      }
    }
    else{
      if(penalized){
        METHOD <- "Penalized weighted Li-Mak test on partial autocorrelations (Gamma Approximation)"
        cor <- Penalized.cor.gof(x, lag=lag, fitdf=fitdf, type="WLM.pacf", alpha=alpha)
        obs <- cor$r.pen
      }
      else{
        METHOD <- "Weighted Li-Mak test on partial autocorrelations (Gamma Approximation)"
        cor <- stats::pacf(x, lag.max = lag, plot=FALSE, na.action=na.pass)
        obs <- cor$acf[1:lag][k]
      }
    }

    weights <- (lag - (fitdf + 1):lag + (fitdf + 1))/lag
    STATISTIC <- n * sum(weights * obs^2)
    if (penalized){
      names(STATISTIC) <- "Penalized weighted X-squared on Squared Residuals for fitted ARCH process"
    }
    else{
      names(STATISTIC) <- "Weighted X-squared on Squared Residuals for fitted ARCH process"
    }

    shape <- (3/4) * (lag + fitdf + 1)^2 *
      (lag - fitdf)/(2 * lag^2 + 3 * lag
                     + 2 * lag * fitdf + 2 * fitdf^2 + 3 * fitdf + 1)
    scale <- (2/3) * (2 * lag^2 + 3 * lag + 2 * lag * fitdf +
                        2 * fitdf^2 + 3 * fitdf + 1)/(lag * (lag + fitdf + 1))
    PARAMETER <- c(Shape = shape, Scale = scale)
    PVAL <- 1 - pgamma(STATISTIC, shape = shape, scale = scale)
  }
  else {
    if (type == "correlation"){
      if(penalized){
        METHOD <- "Penalized Li-Mak test on autocorrelations (Chi-Squared Approximation)"
        cor <- Penalized.cor.gof(x, lag=lag, fitdf=fitdf, type="LM", alpha=alpha)
        obs <- cor$r.pen
      }
      else{
        METHOD <- "Li-Mak test on autocorrelations (Chi-Squared Approximation)"
        cor <- stats::acf(x, lag.max = lag, type="correlation", plot=FALSE, na.action=na.pass)
        obs <- cor$acf[2:(lag + 1)][k]
      }
    }
    else{
      if(penalized){
        METHOD <- "Penalized Li-Mak test on partial autocorrelations (Chi-Squared Approximation)"
        cor <- Penalized.cor.gof(x, lag=lag, fitdf=fitdf, type="LM.pacf", alpha=alpha)
        obs <- cor$r.pen
      }
      else{
        METHOD <- "Li-Mak test on partial autocorrelations (Chi-Squared Approximation)"
        cor <- stats::pacf(x, lag.max = lag, plot=FALSE, na.action=na.pass)
        obs <- cor$acf[1:lag][k]
      }
    }

    weights <- rep(1, (lag - fitdf))
    STATISTIC <- n * sum(weights * obs^2)
    if (penalized){
      names(STATISTIC) <- "Penalized X-squared on Squared Residuals for fitted ARCH process"
    }
    else{
      names(STATISTIC) <- "X-squared on Squared Residuals for fitted ARCH process"
    }
    PARAMETER <- c(df = lag - fitdf)
    PVAL <- 1 - pchisq(STATISTIC, lag - fitdf)
  }


  names(PVAL) <- "Approximate p-value"
  structure(list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = PVAL,
                 method = METHOD,
                 data.name = DNAME),
            class = "htest")
}


# getpvalue <- function(ALPHA,dat,LAG,TYPE,FITDF,SQUARED=FALSE,LOG.SQUARED=FALSE, ABSOLUTE=FALSE,WEIGHTED=TRUE){	
#   return(as.numeric(Penalized.Box.test(x=dat,lag=LAG,type=TYPE,fitdf=FITDF,squared=SQUARED,log.squared=LOG.SQUARED, absolute=ABSOLUTE,weighted=WEIGHTED,alpha=ALPHA)$p.value))
# }
# 
# Penalized.Box.test.pvalue <- function(x, lag = 1, type = c("Ljung-Box", "Box-Pierce", "Monti"),
#                                         fitdf = 0, squared = FALSE, log.squared = FALSE, absolute = FALSE,
#                                         weighted = TRUE, penalized = TRUE) {
# alpha <- 0.05
# if(penalized){
#   pnot <- as.numeric(Penalized.Box.test(x,lag,type,fitdf,squared,log.squared,absolute,weighted,penalized=FALSE)$p.value)
# if(pnot<1e-12) 
# return(Penalized.Box.test(x,lag,type,fitdf,squared,log.squared,absolute,weighted,penalized=FALSE))
# if(pnot>.95)
# return(Penalized.Box.test(x,lag,type,fitdf,squared,log.squared,absolute,weighted,penalized=FALSE))
# low <- pnot/4
# upp <- min(4*pnot,.999)
# step <- (upp-low)/1000
# alphas <- seq(low,upp,step)
# tmp <- sapply(alphas,getpvalue,dat=x,LAG=lag,TYPE=type,FITDF=fitdf,SQUARED=squared,LOG.SQUARED=log.squared,ABSOLUTE=absolute,WEIGHTED=weighted)
#                 	index=order(abs(tmp-alphas))[1]
#                 	alpha=max(tmp[index],min(alphas))
# }
# return(Penalized.Box.test(x, lag, type,
#                           fitdf , squared, log.squared, absolute,
#                           weighted, penalized,alpha))}
