######
#' @title Penalized.LM.test
#'
#' @description The Penalized Li-Mak type test
#'
#' @param x The residuals or initial data which is a numeric vector or univariate time series
#' @param  h.t The sample conditional variances
#' @param lag The lag which the statistic will be based on to calculate penalized acf/pacf, default = 10
#' @param type  The type of test based on acf or pacf, default="correlation"
#' @param fitdf The number of degrees of freedom to be subtracted if x is a series of residuals, default = 0
#' @param alpha The nominal level selected, default = 0.05
#' @param penalized  If FALSE, perform the unpenalized tests using the sample ACF/PACF estimates, default=TRUE
#' @param weighted   -- If FALSE, perform the Li-Mak test, default=TRUE.
#'
#' @return A list with class "htest" for penalized/original Li-Mak tests containing
#'    statistics, parameters, and p value.
#'
#' @examples
#' \dontrun{
#' library(fGarch)
#' spec <- garchSpec(model=list(ar=0.2,omega=0.2, alpha=c(0.2), beta=0))
#' x <- garchSim(spec, n=100)
#' tmp <- garchFit(~arma(1,0)+garch(1,0), data=x, trace=FALSE)
#' res <- residuals(tmp)
#' h.t <- attr(tmp, "h.t")
#' out_test <- Penalized.LM.test(x = res, h.t=h.t, lag=10, type="correlation", fitdf=1, penalized=TRUE)
#' }
#' @export
#####

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

