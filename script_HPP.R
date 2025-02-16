#########################################################################################################
#This file contains the major R code for the test plans based on HPP model in
#paper "Planning Reliability Assurance Tests for Autonomous Vehicles Based on Disengagement Events Data"
#by Simin Zheng, Lu Lu, Yili Hong, and Jian Liu
#########################################################################################################
setwd("/your directory path")
library(plotrix)
source("sharedfuns.r")

#################################
###Apply the functions for HPP###
#################################
source("func_HPP.R")
## Explore a set of test plans for a range of test duration (days) and maximum number of allowable failures 
load("waymo.bayes.we.altp.gibbs.obj")
dmv.hpp.df=tp.eval.hpp(ths=seq(200,3650,by=1),
                       cs=seq(0,50,by=5),
                       m1=0.016,m0=0.013,
                       pars.drs=waymo.bayes.we.altp.gibbs.obj$post.mat,
                       mpd.test=0.21,mpd.use=0.21)

## Produce test plans for HPP
## Need to create a new folder called "figures" under the current working directory
# CR vs PR (Top Left)
fig.paper(filename="HPP_figure_one_CR_PR_total_test_dur",type="pdf")
plot.CR.PR(critvals=dmv.hpp.df,
           ths=seq(200,3650,by=1),
           cs=seq(0,50,by=5),
           x_start=0,x_end=0.15,y_start=0,y_end=0.35,
           t=c(500,1000,1500,2500,3500),
           pt=c(0,15,19,17,18),
           legend_labels=c(expression(paste(tau, "=500")),expression(paste(tau, "=1000")),
                           expression(paste(tau, "=2000")),expression(paste(tau, "=3000")),
                           expression(paste(tau, "=4000"))),
           col_labels= c("c=0","c=10","c=20","c=30","c=40","c=50"),
           cl_xleft=0,cl_ybottom=0.05,cl_xright=0.0015,cl_ytop=0.10)
dev.off()

# TH vs AP (Top Right)
fig.paper(filename="HPP_figure_two_TH_AP_total_test_dur",type="pdf")
plot.TH.AP(critvals=dmv.hpp.df,
           ths=seq(200,3650,by=1),
           cs=seq(0,50,by=5),
           x_adjust=530,
           cr=c(0.04,0.06,0.08,0.10,0.12),pt1=c(1:5),
           pr=c(0.10,0.16,0.22,0.28,0.34),pt2=c(17,18,20,15,0),
           l_x=3050,l_y=1.00,
           legend_labels=c('CR=0.04','CR=0.06','CR=0.08','CR=0.10','CR=0.12','PR=0.10','PR=0.16','PR=0.22','PR=0.28','PR=0.34'),
           width=400,
           col_labels= c("c=0","c=10","c=20","c=30","c=40","c=50"),
           cl_xleft=3720.97,cl_ybottom=0.67,cl_xright=3760,cl_ytop=0.82)
dev.off()

# TH vs CR (Bottom Left)
fig.paper(filename="HPP_figure_three_TH_CR_total_test_dur",type="pdf")
plot.TH.CR(critvals=dmv.hpp.df,
           ths=seq(200,3650,by=1),
           cs=seq(0,50,by=5),
           ystart=0,yend=0.15,
           pr=c(0.10,0.16,0.22,0.28,0.34),
           pt=c(18,17,20,15,0),
           l_x=3250,l_y=0.15,
           legend_labels=c('PR=0.10','PR=0.16','PR=0.22','PR=0.28','PR=0.34'),
           cl_xleft=3300.97,cl_ybottom=0.096,cl_xright=3340,cl_ytop=0.123)
dev.off()

# TH vs PR (Bottom Right)
fig.paper(filename="HPP_figure_four_TH_PR_total_test_dur",type="pdf")
plot.TH.PR(critvals=dmv.hpp.df,
           ths=seq(200,3650,by=1),
           cs=seq(0,50,by=5),
           ystart=0,yend=0.36,
           cr=c(0.04,0.06,0.08,0.10,0.12),
           pt=c(1:5),
           x_adjust=500,l_x=3700,l_y=0.12,
           legend_labels=c('CR=0.04','CR=0.06','CR=0.08','CR=0.10','CR=0.12'),
           cl_xleft=3800.97,cl_ybottom=0.0,cl_xright=3840,cl_ytop=0.057,
           col_labels= c("c=0","c=10","c=20","c=30","c=40","c=50"))
dev.off()

##########################################
###Apply Pareto Front functions for HPP###
##########################################
dmv.hpp.df_each_c=tp.eval.hpp(ths=seq(200,3650,by=1),
                              cs=seq(0,50,by=1),
                              m1=0.016,m0=0.013,
                              pars.drs=waymo.bayes.we.altp.gibbs.obj$post.mat,
                              mpd.test=0.21,mpd.use=0.21)

# Pareto Front optimization
source("func_Pareto_Front.R")
tp.hpp.pf.df_check=tp.dat(dmv.hpp.df_each_c,ths=seq(200,3650,by=1),cs=seq(0,50,by=1))

# Find the Pareto optimal test plans for a selected cut off based on consumer's risk
tp.hpp.pf_check=pf.cond.cr(tp.set=tp.hpp.pf.df_check,cr.cut=0.086)
tp.hpp.pf_check[[2]] #computation time for the HPP model

# Plot the Pareto optimal test plans for the HPP model
fig.paper(filename="HPP_Pareto_Front_for_each_c",type="pdf")
plot.pf(tp.hpp.pf_check[[1]])
dev.off()