# This code was used to conduct a simulation study to compare the performance
# of the estimation method for different choices of the tuning parameter (grid length).

#load function to simulate survival data
source("sim.survdata.R")

# install and load package CPsurv from GitHub
devtools::install_github("stefkruegel/CPsurv")
library(CPsurv)

# load further packages
library(foreach)
library(doParallel)
library(snow)


# fixed parameters:
nsamp <- 1000   # number of simulation samples
shape <- 0.44
scale <- 100
cpmax <- 360
censpoint <- 540
B.correct = 49
B = 999
cores <- min(19, parallel::detectCores())   #CPU-cores

# varying parametes for simulating data:
n <- c(1000, 5000)   # number of observations
cens <- c("random", "type1", "no")   # censoring type
jump <- c(1, 0)   #jump in hazard rate at changepoint
chp <- c(55, 90)  # changepoint

# all possible parameter combinations:
param <- as.data.frame(expand.grid(n, chp, jump, cens))
names(param) <- c("n", "chp", "jump", "cens")
param$cens <- as.character(param$cens)

# parameter for estimation: interval width
intwd <- c(3, 5, 10, 20)


set.seed(20160701)
seeds <- sample.int(1e8, size = nsamp)


results <- list()

for(i in 1:nrow(param)){

  cl <- makeCluster(cores)
  registerDoParallel(cl)
  resultlist <- foreach(idx = 1:nsamp, .packages = "CPsurv") %dopar% {
    set.seed(seeds[idx])

    # simulate main dataset of length 'nobs' with parameters param[i,]
    simdata <- sim.survdata(nobs = param$n[i], changeP = param$chp[i],
                            shape = shape, scale = scale, jump = param$jump[i],
                            censoring = param$cens[i], censpoint = censpoint)

    #---------------------------------------------------------------------------

    resmat <- matrix(nrow = length(intwd), ncol = 10)
    for(j in 1:length(intwd)){

      est <- cpsurv(time = simdata$time, event = simdata$event,
                    cpmax = cpmax, intwd = intwd[j],
                    censoring = param$cens[i], censpoint = censpoint,
                    biascorrect = FALSE, parametric = FALSE,
                    B.correct = B.correct, boot.ci = TRUE, B = B,
                    parallel = FALSE)
      
      tab <- table(cut(simdata$time[simdata$event == 1], est$lower.lim))
      tab2 <- table(cut(simdata$time[simdata$event == 1], 
                        est$lower.lim[est$lower.lim <= param$chp[i]]))

      cilow <-    est$ci.normal[1]
      ciup <-     est$ci.normal[2]
      perclow <-  est$ci.percent[1]
      percup <-   est$ci.percent[2]
      std <-      est$sd
      min.events <- min(tab)
      mean.events <- mean(tab)
      min.events.chp <- min(tab2)

      resmat[j, ] <- c(est$cp, est$pv.mean, cilow, ciup, perclow, percup, std,
                       min.events, mean.events, min.events.chp)
    }
    rownames(resmat) <- intwd  
    resmat
  }
  stopCluster(cl)

  #============================================================================
  output <- matrix(nrow = 0, ncol = 10)
  colnames(output) <- c("cp", "pv.mean", "ci1", "ci2", "perc1", "perc2", "sd", 
                        "min.events", "mean.events", "min.events.chp")
  for(m in 1:length(resultlist)){
     output <- rbind(output, resultlist[[m]])
  }

  names <- paste("cens=", param$cens[i], "; jump=", param$jump[i],
                 "; chp=", param$chp[i], "; n=", param$n[i], sep="")
  results[[names]] <- output
  counter <- c(i)
  save(counter, file="counter.rda")
  save(results, file="simulation_results.rda")
}



##############################################################################
### processing the results
##############################################################################

# function for root mean squared error (RMSE)
rmse <- function(x, i){
  if(i %in% c(1,2,5,6,9,10,13,14,17,18,21,22)){
    tau <- 55
  }else{
    tau <- 90
  }
  ret <- round(sqrt( mean((x-tau)^2, na.rm = TRUE) ),2)
  ret
}

# function for mean average distance (MAD)
mad <- function(x, i){
  if(i %in% c(1,2,5,6,9,10,13,14,17,18,21,22)){
    tau <- 55
  }else{
    tau <- 90
  }
  ret <- round( mean(abs(x-tau), na.rm = TRUE) ,2)
  ret
}

# function for calculating the coverage of confidence or percentile intervals
cover <- function(lower, upper, i){
  if(i %in% c(1,2,5,6,9,10,13,14,17,18,21,22)){
    tau <- 55
  }else{
    tau <- 90
  }
  abscover <- sum(sapply(seq_along(lower), 
                         function(i) tau >= lower[i] && tau <= upper[i]))
  round(abscover / length(lower), 3)
}


# calculate lengths of confidence and percentile intervals
for(i in 1:24){
  results[[i]] <- cbind(results[[i]], results[[i]][, "ci2"] - results[[i]][, "ci1"])
  results[[i]] <- cbind(results[[i]], results[[i]][, "perc2"] - results[[i]][, "perc1"])
  colnames(results[[i]])[11:12] <- c("ci.len", "perc.len")
}

#===========================================================================
names.est <- rownames(results$`cens=random; jump=1; chp=55; n=1000`[1:24,])

measures <- matrix(NA, nrow = (length(results) * 4), ncol = 11)
colnames(measures) <- c("median", "mean", "mad", "rmse", "ci.length", "ci.cover", 
                        "perc.length", "perc.cover", "min.events", "mean.events",
                        "min.events.beforecp")

for(i in (seq_along(results)-1)){
  measures[(1:4)+i*4, 1] <- sapply(1:4, function(x) round(median(results[[i+1]][x+((1:1000)-1)*4, "cp"], na.rm = TRUE)))
  measures[(1:4)+i*4, 2] <- sapply(1:4, function(x) round(mean(results[[i+1]][x+((1:1000)-1)*4, "cp"], na.rm = TRUE)))
  measures[(1:4)+i*4, 3] <- sapply(1:4, function(x) mad(results[[i+1]][x+((1:1000)-1)*4, "cp"], i+1))
  measures[(1:4)+i*4, 4] <- sapply(1:4, function(x) rmse(results[[i+1]][x+((1:1000)-1)*4, "cp"], i+1))
  measures[(1:4)+i*4, 5] <- sapply(1:4, function(x) round(mean(results[[i+1]][x+((1:1000)-1)*4, "ci.len"], na.rm = TRUE)))
  measures[(1:4)+i*4, 6] <- sapply(1:4, function(x) cover(results[[i+1]][x+((1:1000)-1)*4, "ci1"], 
                                                          results[[i+1]][x+((1:1000)-1)*4, "ci2"], i+1))
  measures[(1:4)+i*4, 7] <- sapply(1:4, function(x) round(mean(results[[i+1]][x+((1:1000)-1)*4, "perc.len"], na.rm = TRUE)))
  measures[(1:4)+i*4, 8] <- sapply(1:4, function(x) cover(results[[i+1]][x+((1:1000)-1)*4, "perc1"], 
                                                          results[[i+1]][x+((1:1000)-1)*4, "perc2"], i+1))
  measures[(1:4)+i*4, 9] <- sapply(1:4, function(x) round(mean(results[[i+1]][x+((1:1000)-1)*4, "min.events"], na.rm = TRUE)))
  measures[(1:4)+i*4, 10] <- sapply(1:4, function(x) round(mean(results[[i+1]][x+((1:1000)-1)*4, "mean.events"], na.rm = TRUE)))
  measures[(1:4)+i*4, 11] <- sapply(1:4, function(x) round(mean(results[[i+1]][x+((1:1000)-1)*4, "min.events.chp"], na.rm = TRUE)))
}

#===========================================================================

# convert matrix to LaTeX table

library(xtable)
measures1 <- data.frame(row=rep(names.est, 4), measures)
xmeasures <- xtable(measures1)

addtorow <- list()
addtorow$pos <- list()
addtorow$command <- vector("numeric")

for(i in 1:24){
  names <- names(results)[i]
  addtorow$pos[[i]] <- c(4*i-4)
  addtorow$command[i] <- paste0("\\newpage \n \\multicolumn{12}{l}{\\textit{", names, "}}\\\\ \n")
}

print(xmeasures, add.to.row = addtorow, include.rownames = FALSE,
      include.colnames = TRUE, tabular.environment="longtable", floating=FALSE)

