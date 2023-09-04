######
#' @title Penalized.Box.test
#'
#' @description The Penalized (Weighted) Box-Pierce, (Weighted) Ljung-Box, (Weighted) Monti, Weighted McLeod-Li type test for fitted ARMA and detection of nonlinear Models
#'
#' @param x The residuals or initial data which is a numeric vector or univariate time series
#' @param lag The lag which the statistic will be based on to calculate penalized acf/pacf, default = 10
#' @param type the type of test, default = "Ljung-Box"
#' @param fitdf The number of degrees of freedom to be subtracted if x is a series of residuals, default = 0
#' @param alpha The nominal level selected, default = 0.05
#' @param squared Take the squared value of residuals to detect nonlinear processes, default=FALSE
#' @param log.squared Take the log of the squared residuals to detect nonlinear processes, default=FALSE
#' @param absolute Take the absolute value of the residuals to detect nonlinear processes, default=FALSE
#' @param penalized  If FALSE, perform the unpenalized tests using the sample ACF/PACF estimates, default=TRUE
#' @param weighted   If FALSE, perform the original Box-Pierce, Ljung-Box or Monti or McLeod Li type if using a transformation, default=TRUE.
#'
#' @return A list with class "htest" for penalized/traditional (Weighted) Box-Pierce,
#'      (Weighted) Ljung-Box, (Weighted) Monti tests
#'      containing statistics, parameters, and p value.
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.9))
#' ## Fit AR(1) model
#' p.fit <- 1
#' q.fit <- 0
#' res <- arima(data, order=c(p.fit, 0, q.fit))$resid
#' out_test1 <- Penalized.Box.test(x = res, lag=10, fitdf=1, type="Ljung-Box", weighted=TRUE, penalized=TRUE)
#'
#' ## Fit ARMA(1,1) model
#' p.fit <- 1
#' q.fit <- 1
#' res <- arima(data, order=c(p.fit, 0, q.fit))$resid
#' out_test2 <- Penalized.Box.test(x = res, lag=10, fitdf=1, type="Ljung-Box", weighted=TRUE, penalized=TRUE)
#' }
#' @export
#####

Penalized.Box.test <- function(x, lag = 1, type = c("Ljung-Box", "Box-Pierce", "Monti"),
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
