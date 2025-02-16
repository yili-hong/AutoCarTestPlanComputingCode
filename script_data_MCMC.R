#########################################################################################################
#This file contains the major R code for generating the Waymo MCMC in
#paper "Planning Reliability Assurance Tests for Autonomous Vehicles Based on Disengagement Events Data"
#by Simin Zheng, Lu Lu, Yili Hong, and Jian Liu
#########################################################################################################
setwd("/your directory path")

source("sharedfuns.r")
library(Rcpp)
sourceCpp("functions.cpp")

## Read in the data from csv files
waymo.dat<-list(
  event.dat=read.csv("data/event.dat.csv"),
  mile.dat=read.csv("data/mile.dat.csv"),
  initial.time=as.Date("2017-11-30"),
  time.int.tab=read.csv("data/time.int.tab.csv")
  
)
waymo.dat$event.dat$Date<-as.Date(waymo.dat$event.dat$Date, format = "%m/%d/%y")
waymo.dat$time.int.tab$start.time<-as.Date(waymo.dat$time.int.tab$start.time, format = "%Y-%m-%d")
waymo.dat$time.int.tab$end.time<-as.Date(waymo.dat$time.int.tab$end.time, format = "%Y-%m-%d")

waymo.dat.obj=data.setup(dat=waymo.dat, manufacturer="Waymo")

## Try Weibull alternative parametrization
# fit MLE 
minus.loglik.weibull.altp(dat=waymo.dat.obj$count.obj, 
                          pars=c(4, 5, 0))

fit.we.altp.waymo=lifetime.mle(dat=waymo.dat.obj$count.obj, 
                               minusloglik=minus.loglik.weibull.altp, 
                               starts=c(4, 5, 0), 
                               method = "Nelder-Mead", 
                               control=list(trace=T, maxit=2000))

## Run Bayesian with normal priors.
# waymo MCMC
waymo.bayes.we.altp.gibbs.obj=Bayes.wei.altp.fit(dat=waymo.dat.obj$count.obj, 
                                                 init=fit.we.altp.waymo$coef, 
                                                 B=1.1e6, 
                                                 B.sample=.1e6, 
                                                 b.sample=1e3, 
                                                 prior.inf=c(4.6, 3, 6.7, 3, 0, 3))

# save the above Gibbs samples to your directory for further analysis
save(waymo.bayes.we.altp.gibbs.obj, file="waymo.bayes.we.altp.gibbs.obj")