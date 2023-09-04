######
#' @title Penalized.cor.gof
#'
#' @description The Penalized ACF/PACF estimator for time series goodness of fit
#'
#' @param x A numeric vector or univariate time series
#' @param lag The lag which the statistic will be based on to calculate penalized acf/pacf, default = 10
#' @param type The type of test which penalized acf/pacf to be used, default = "LB"
#' @param fitdf The number of degrees of freedom to be subtracted if x is a series of residuals, default = 0
#' @param alpha The nominal level selected, default = 0.05
#'
#' @return An object of class penalized acf/pacf estimation for time series goodness of fit tests
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.9))
#' pen_cor <- Penalized.cor.gof(x = data, lag=10, fitdf=1, type="LB", alpha=0.05)
#' }
#' @export
#####
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


