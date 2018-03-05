### R-Code to reproduce simulation results and data application in "Nonparametric change point estimation for survival distributions with a partially constant hazard rate"

The paper mentioned above, which contains the underlying methods of ```CPsurv```, will be published in the journal "Lifetime Data Analysis" next time. It contains an extensive simulation study which was partly conducted with an unpublished version of ```CPsurv``` where different variations of the estimation method were included. The last part of the simulation study, where several values of the tuning parameter (grid length) were compared, is feasible with the current version of the package. 

Run the ```simulation-study.R``` script to reproduce the results published in section 4.4 of the paper. Please note the extended runtime for estimation with bias correction, particularly in combination with calculation of confidence intervals (because variance estimation of the changepoint is based on a bootstrap procedure). 

Run the ```data-examples.R``` script to reproduce the estimation results in section 5 of the paper. The three data sets include real survival times of ICU (intensive care unit) patients. 

The function ```sim.survdata``` is needed to simulate survival data for the simulation study.
