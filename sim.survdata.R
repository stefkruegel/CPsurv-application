#' Simulate Survival Data with Changepoint
#'
#' Simulates Weibull distributed survival data with changepoint above which
#' hazard rate is constant.
#'
#' @inheritParams cpest
#' @param changeP Changepoint.
#' @param shape Shape parameter of Weibull distribution.
#' @param scale Scale parameter of Weibull distribution.
#' @param censoring Logical; if \code{TRUE}, censored data are generated.
#' @param censpoint Censoring point for Type I censoring.
#' @param times.int Logical; if \code{TRUE}, returned survival times are
#'   integers
#' @param parametric Logical; if \code{TRUE}, survival times are generated
#'   parametrically by inverse transform sampling; otherwise Kaplan-Meier is
#'   used for simulation

#' @return A dataset with survival times and corresponding censoring status
#'   ('event').
sim.survdata <-function(nobs, changeP, shape = 0.44, scale = 100, jump = FALSE,
                        censoring = "random", censpoint = 540, times.int = FALSE,
                        rate.cens = 0.0025){

  # rate of exponential distr.
  rateE <- shape / scale * (changeP / scale) ^ (shape - 1)
  if(jump) rateE <- rateE / 2

  # Survivorfct. at changepoint
  St <- pweibull(q = changeP, shape = shape, scale = scale, lower.tail = FALSE)

  u <- runif(nobs)
  times.sim <- numeric(nobs)
  times.sim[u > St] <- qweibull(p = u[u > St], shape = shape, scale = scale,
                                lower.tail = FALSE)
  times.sim[u <= St] <- changeP + qexp(p = u[u <= St] /
                                         pweibull(q = changeP, shape = shape,
                                                  scale = scale,
                                                  lower.tail = FALSE),
                                       rate = rateE, lower.tail = FALSE)


  data <- data.frame(time = times.sim, event = rep(1, nobs))

  #-----------------censoring--------------------
  if(censoring == "random"){
    censtimes <- rexp(nobs, rate.cens)
    data$event <- as.numeric(data$time <= censtimes)
    data$time <- pmin(data$time, censtimes)
  }else if(censoring == "type1"){
    data$event[data$time >= censpoint] <- 0
    data$time[data$time >= censpoint] <- censpoint
  }

  if(times.int){
    data$time <- ceiling(data$time)
  }
  data
}
