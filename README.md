survPCA: A package for estimating martingale and counting process covariance matrix at a specified time
=======================================================================================================

The `survPCA` package estimates the covariance matrix for martingales and counting processes both in the semi-competing risk setting and the no competing risk setting.

Installation
------------

``` install
devtools::install_github("blangworthy/survPCA")
library(survPCA)
```

Examples
--------

``` r
library(survPCA)
library(MASS)

#####Example without competing risks
###In this example we simulate 4 non-competing events and calculate the martingale and counting process covariance and correlation matrices. All events have exponential 1 distribution and we estimate the covariance/correlation matrices at timepoint 1.

n <- 200 # sample size
sigma <- toeplitz(c(1,0.6,0.4,0.2)) # Correlation structure for gaussian copula
samps <- mvrnorm(n,rep(0,4),sigma) #Simulate from multivariate normal
psamps <- pnorm(samps) #Transform marginals to uniform
expsamps <- -log(1-psamps) #Transform marginals to exponential(1), this is failure time


cens <- runif(n,min=0,max=4) #Censoring time

eventsamps <-apply(expsamps,2,function(x)ifelse(x< cens,1,0)) #Censoring indicator

observed_time <- apply(expsamps,2,function(x) pmin(x,cens)) #Minimum of failur time and censoring time

cov_est_data <- cbind(observed_time,eventsamps) #Data in format for CovMatNoComp function

t <- 1 #Time to estimate martingale/counting process covariance and correlation matrices

cov_out <- CovMatNoComp(data = cov_est_data,p = 4,t = t) #Gives output of martingale/counting process correlation and covariance matrices
pc_mart_cor <- eigen(cov_out$CorMart)$vectors # PC directions for martingale correlation matrix

######Example with competing risks
###In this set up there are 4 non-competing events which we will calculate the martingale/counting process correlation and covariance matrices. 
###In addition there is one competing event as well as a censoring time. We will use the same sample size and censoring time as previous example.

sigma <- rbind(cbind(toeplitz(c(1,0.6,0.4,0.2)),rep(0.1,4)),c(rep(0.1,4),1))#This has one extra column/row representing the competing risk time
samps_comp <- mvrnorm(n,rep(0,5),sigma) #Simulate from multivariate normal
psamps_comp <- pnorm(samps_comp) #Transform marginals to uniform
expsamps_comp <- -log(1-psamps_comp) #Transform marginals to exponential(1), this is failure time


noncompsamps <- expsamps_comp[,(1:4)]#Failure ime for non-competing events
compsamps <- expsamps_comp[,5]#Failure time for competing event

causeind <- apply(noncompsamps,2,function(x) ifelse(x<compsamps,1,2))#Cause indicator, 1 if non-competing event, 2 if competing event
failtime <- apply(noncompsamps,2,function(x) pmin(x,compsamps))#Minimum of non competing event and competing event

censind <- apply(failtime,2,function(x) ifelse(x<cens,1,0))#Indicator of whether censored

observed_time <- apply(failtime,2,function(x) pmin(x,cens)) #Minimum of failure time and censoring time

eventsamps <- causeind*censind #0 if censored, 1 if failure time observed and non-competing event, 2 if failure time observed and competing event

cov_est_data_comp <- cbind(observed_time,eventsamps) #Data in format for CovMatComp function
cov_out_comp <- CovMatComp(data = cov_est_data_comp,causes = rep(1,4),p = 4, t=t) #Gives output of martingale/counting process correlation and covariance matrices
pc_mart_cor_comp<- eigen(cov_out_comp$CorMart)$vectors # PC directions for martingale correlation matrix
```
