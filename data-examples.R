# install package CPsurv from GitHub
devtools::install_github("stefkruegel/CPsurv")

library(CPsurv)

################################################################################
### Data set 1: Critically ill patients
################################################################################

load("dataset1.RData")

data1result <- list()

### estimate changepoint (raw estimate)

# 10 days grid length
data1result$raw10 <- cpsurv(dataset1$time, dataset1$event, intwd = 10, cpmax = 300,
                         censoring = "type1", censpoint = 730,
                         boot.ci = TRUE, seed = 20160501, cores = 6)
summary(data1result$raw10)

# 20 days grid length
data1result$raw20 <- cpsurv(dataset1$time, dataset1$event, intwd = 20, cpmax = 300,
                         censoring = "type1", censpoint = 730,
                         boot.ci = TRUE, seed = 20160501, cores = 6)
summary(data1result$raw20)

### estimate changepoint with nonparametric median bias correction
# (note extended runtime because of bias correction: 
# about 30 minutes using 6 cores, Intel i7 processor)

# 10 days grid length
data1result$mbc10 <- cpsurv(dataset1$time, dataset1$event, intwd = 10, cpmax = 300, 
                         censoring = "type1", censpoint = 730, biascorrect = TRUE,
                         boot.ci = TRUE, opt.start = c(0.8, 2500), seed = 20160501, 
                         cores = 6)
summary(data1result$mbc10)

# 20 days grid length
data1result$mbc20 <- cpsurv(dataset1$time, dataset1$event, intwd = 20, cpmax = 300, 
                         censoring = "type1", censpoint = 730, biascorrect = TRUE,
                         boot.ci = TRUE, opt.start = c(0.8, 2500), seed = 20160501, 
                         cores = 6)
summary(data1result$mbc20)

# estimate of beta (expectation of p-values for intervals above the changepoint)
data1result$raw10$pv.mean
data1result$raw20$pv.mean


################################################################################
### Data set 2: Patients after severe trauma injury
################################################################################

load("dataset2.RData")

data2result <- list()

### estimate changepoint (raw estimate)

# 10 days grid length
data2result$raw10 <- cpsurv(dataset2$time, dataset2$event, intwd = 10, cpmax = 700,
                            boot.ci = TRUE, seed = 20160501, cores = 6)
summary(data2result$raw10)

# 20 days grid length
data2result$raw20 <- cpsurv(dataset2$time, dataset2$event, intwd = 20, cpmax = 700,
                            boot.ci = TRUE, seed = 20160501, cores = 6)
summary(data2result$raw20)


### estimate changepoint with nonparametric median bias correction

# 10 days grid length
data2result$mbc10 <- cpsurv(dataset2$time, dataset2$event, intwd = 10, cpmax = 700, 
                            biascorrect = TRUE, boot.ci = TRUE, 
                            opt.start = c(0.2, 200000000), seed = 20160501, 
                            cores = 6)
summary(data2result$mbc10)

# 20 days grid length
data2result$mbc20 <- cpsurv(dataset2$time, dataset2$event, intwd = 20, cpmax = 700, 
                            biascorrect = TRUE, boot.ci = TRUE, 
                            opt.start = c(0.2, 200000000), seed = 20160501, 
                            cores = 6)
summary(data2result$mbc20)


# estimate of beta (expectation of p-values for intervals above the changepoint)
data2result$raw10$pv.mean
data2result$raw20$pv.mean



################################################################################
### Data set 3: Postoperative survival of patients after a partial hepatectomy
################################################################################

load("dataset3.RData")

data3result <- list()

### estimate changepoint (raw estimate)

# 10 days grid length
data3result$raw10 <- cpsurv(dataset3$time, dataset3$event, intwd = 10, cpmax = 200,
                            boot.ci = TRUE, seed = 20160501, cores = 6)
summary(data3result$raw10)

# 20 days grid length
data3result$raw20 <- cpsurv(dataset3$time, dataset3$event, intwd = 20, cpmax = 200,
                            boot.ci = TRUE, seed = 20160501, cores = 6)
summary(data3result$raw20)


### estimate changepoint with nonparametric median bias correction

# 10 days grid length
data3result$mbc10 <- cpsurv(dataset3$time, dataset3$event, intwd = 10, cpmax = 200, 
                            biascorrect = TRUE, boot.ci = TRUE, 
                            opt.start = c(0.8, 2500), seed = 20160501, cores = 6)
summary(data3result$mbc10)

# 20 days grid length
data3result$mbc20 <- cpsurv(dataset3$time, dataset3$event, intwd = 20, cpmax = 200, 
                            biascorrect = TRUE, boot.ci = TRUE, 
                            opt.start = c(0.8, 2500), seed = 20160501, cores = 6)
summary(data3result$mbc20)


# estimate of beta (expectation of p-values for intervals above the changepoint)
data3result$raw10$pv.mean
data3result$raw20$pv.mean



