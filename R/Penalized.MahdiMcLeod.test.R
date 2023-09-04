######
#' @title Penalized.MahdiMcLeod.test
#'
#' @description The Penalized MahdiMcLeod test for fitted ARMA and detection of nonlinear Models
#'
#' @param x The residuals or initial data which is a numeric vector or univariate time series
#' @param lag The lag which the statistic will be based on to calculate penalized acf/pacf, default = 10
#' @param fitdf The number of degrees of freedom to be subtracted if x is a series of residuals, default = 0
#' @param alpha The nominal level selected, default = 0.05
#' @param squared Take the squared value of residuals to detect nonlinear processes, default=FALSE
#' @param log.squared Take the log of the squared residuals to detect nonlinear processes, default=FALSE
#' @param absolute Take the absolute value of the residuals to detect nonlinear processes, default=FALSE
#' @param penalized  If FALSE, perform the unpenalized tests using the sample ACF/PACF estimates, default=TRUE
#'
#' @return A list with class "htest" for penalized/tranditional MahdiMcLeod test results
#'     containing statistics, parameters, and p value.
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.9))
#' ## Fit AR(1) model
#' p.fit <- 1
#' q.fit <- 0
#' res <- arima(data, order=c(p.fit, 0, q.fit))$resid
#' out_test1 <- Penalized.MahdiMcLeod.test(x = res, lag=10, fitdf=1, penalized=TRUE)
#' #'
#' ## Fit ARMA(1,1) model
#' p.fit <- 1
#' q.fit <- 1
#' res <- arima(data, order=c(p.fit, 0, q.fit))$resid
#' out_test2 <- Penalized.MahdiMcLeod.test(x = res, lag=10, fitdf=1, penalized=TRUE)
#' }
#' @export
#####

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

