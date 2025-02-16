#########################################################################################################
#This file contains the major R code for the test plans based on NHPP model in
#paper "Planning Reliability Assurance Tests for Autonomous Vehicles Based on Disengagement Events Data"
#by Simin Zheng, Lu Lu, Yili Hong, and Jian Liu
##########################################################################################################
setwd("/your directory path")
library(plotrix)
source("sharedfuns.r")

########################
####Apply NHPP func.####
########################
source("func_NHPP.R")
load("waymo.bayes.we.altp.gibbs.obj")
## Testing period: 1 year
## For single vehicle test
# CR vs PR
dmv.nhpp.df_1.test=tp.eval.nhpp(ths=seq(0.1*365,365,by=5), 
                                cs=seq(0,5,by=1),
                                m1=0.0125,m0=0.009,
                                pars.drs=waymo.bayes.we.altp.gibbs.obj$post.mat,
                                mpd.test=0.20,mpd.use=0.20,dem.dur=2*365, nt=1)
# Visualizations for nt=1
fig.paper(filename="NHPP_figure_one_CR_PR_plot_1",type="pdf")
plot.CR.PR(critvals=dmv.nhpp.df_1.test,
           ths=seq(0.1*365,1*365,by=5),
           cs=seq(0,5,by=1),
           x_start=0,x_end=0.27,y_start=0.1,y_end=0.35,
           t=c(50,100,200,300),pt=c(0,15,19,17,18),
           legend_labels=c(expression(paste(tau[t], "=50")),expression(paste(tau[t], "=100")),
                           expression(paste(tau[t], "=200")),expression(paste(tau[t], "=300"))),
           color_labels=c("c=0","c=1","c=2","c=3","c=4","c=5"),
           cl_xleft=-0.006,cl_ybottom=0.13,cl_xright=-0.001,cl_ytop=0.18)
dev.off()

#(2) TH vs AP
fig.paper(filename="NHPP_figure_two_TH_AP_plot_1",type="pdf")
plot.TH.AP(critvals=dmv.nhpp.df_1.test,
           ths=seq(0.1*365,1*365,by=5),
           cs=seq(0,5,by=1),
           x_adjust=200,
           cr=c(0.15,0.18,0.21,0.24,0.27),pt1=c(1:5),
           pr=c(0.21,0.23,0.25,0.27,0.29),pt2=c(17,18,20,15,0),
           l_x=390,l_y=1.0, 
           legend_label=c('CR=0.15','CR=0.18','CR=0.21','CR=0.24','CR=0.27','PR=0.21','PR=0.23','PR=0.25','PR=0.27','PR=0.29'),
           width=60,color_labels=c("c=0","c=1","c=2","c=3","c=4","c=5"),
           cl_xleft=508,cl_ybottom=0.60,cl_xright=512.3,cl_ytop=0.8)
dev.off()

#(3) TH vs CR
fig.paper(filename="NHPP_figure_three_TH_CR_plot_1",type="pdf")
plot.TH.CR(critvals=dmv.nhpp.df_1.test,
           ths=seq(0.1*365,1*365,by=5),
           cs=seq(0,5,by=1),
           ystart=0,yend=0.28,
           pr=c(0.21,0.23,0.25,0.27,0.29),pt=c(18,17,20,15,0),
           l_x=35,l_y=0.10,
           legend_label=c('PR=0.21','PR=0.23','PR=0.25','PR=0.27','PR=0.29'),
           color_labels=c("c=0","c=1","c=2","c=3","c=4","c=5"),
           cl_xleft=42,cl_ybottom=0,cl_xright=45.8,cl_ytop=0.05)
dev.off()

#(4) TH vs PR
fig.paper(filename="NHPP_figure_four_TH_PR_plot_1",type="pdf")
plot.TH.PR(critvals=dmv.nhpp.df_1.test,
           ths=seq(0.1*365,1*365,by=5),
           cs=seq(0,5,by=1),
           x_adjust=100,ystart=0.1,yend=0.35,
           cr=c(0.15,0.18,0.21,0.24,0.27),pt=c(1:5),
           l_x=390,l_y=0.35,
           legend_label=c('CR=0.15','CR=0.18','CR=0.21','CR=0.24','CR=0.27'),
           color_labels=c("c=0","c=1","c=2","c=3","c=4","c=5"),
           cl_xleft=401.7,cl_ybottom=0.23,cl_xright=406,cl_ytop=0.29)
dev.off()
##########################################################################################################
## For multiple vehicle test
#(1) CR vs PR
dmv.nhpp.df_5=tp.eval.nhpp(ths=seq(0.1*365,1*365,by=5), 
                           cs=seq(0,25,by=1),
                           m1=0.0132,m0=0.009,
                           pars.drs=waymo.bayes.we.altp.gibbs.obj$post.mat,
                           mpd.test=0.2,mpd.use=0.2,dem.dur=2*365, nt=5)

# Visualizations for nt=5
# CR vs PR
fig.paper(filename="NHPP_figure_one_CR_PR_plot_5",type="pdf")
plot.CR.PR(critvals=dmv.nhpp.df_5,
           ths=seq(0.1*365,1*365,by=5),
           cs=seq(0,25,by=1),
           x_start=0,x_end=0.21,y_start=0,y_end=0.38,
           t=c(50,100,200,300),pt=c(0,15,19,17,18),
           legend_labels=c(expression(paste(tau[t], "=50")),expression(paste(tau[t], "=100")),
                           expression(paste(tau[t], "=200")),expression(paste(tau[t], "=300"))),
           color_labels=c("c=0","c=5","c=10","c=15","c=20","c=25"),
           cl_xleft=0,cl_ybottom=0.05,cl_xright=0.003,cl_ytop=0.11)
dev.off()

# TH vs AP
fig.paper(filename="NHPP_figure_two_TH_AP_plot_5",type="pdf")
plot.TH.AP(critvals=dmv.nhpp.df_5,
           ths=seq(0.1*365,1*365,by=5),
           cs=seq(0,25,by=1),
           x_adjust=150,
           cr=c(0.11,0.13,0.15,0.17,0.19),pt1=c(1:5),
           pr=c(0.15,0.20,0.25,0.30,0.35),pt2=c(17,18,20,15,0),
           l_x=380,l_y=1.0, 
           legend_label=c('CR=0.11','CR=0.13','CR=0.15','CR=0.17','CR=0.19','PR=0.15','PR=0.20','PR=0.25','PR=0.30','PR=0.35'),
           width=60,color_labels=c("c=0","c=5","c=10","c=15","c=20","c=25"),
           cl_xleft=408,cl_ybottom=0.60,cl_xright=412.3,cl_ytop=0.80)
dev.off()

# TH vs CR
fig.paper(filename="NHPP_figure_three_TH_CR_plot_5",type="pdf")
plot.TH.CR(critvals=dmv.nhpp.df_5,
           ths=seq(0.1*365,1*365,by=5),
           cs=seq(0,25,by=1),
           ystart=0,yend=0.21,
           pr=c(0.15,0.20,0.25,0.30,0.35),
           pt=c(18,17,20,15,0),
           l_x=300,l_y=0.21,
           legend_label=c('PR=0.15','PR=0.20','PR=0.25','PR=0.30','PR=0.35'),
           color_labels=c("c=0","c=5","c=10","c=15","c=20","c=25"),
           cl_xleft=330,cl_ybottom=0.13,cl_xright=334.3,cl_ytop=0.17)
dev.off()

# TH vs PR
fig.paper(filename="NHPP_figure_four_TH_PR_plot_5",type="pdf")
plot.TH.PR(critvals=dmv.nhpp.df_5,
           ths=seq(0.1*365,1*365,by=5),
           cs=seq(0,25,by=1),
           x_adjust=100,ystart=0,yend=0.36,
           cr=c(0.11,0.13,0.15,0.17,0.19),pt=c(1:5),
           l_x=380,l_y=0.36,
           legend_label=c('CR=0.11','CR=0.13','CR=0.15','CR=0.17','CR=0.19'),
           color_labels=c("c=0","c=5","c=10","c=15","c=20","c=25"),
           cl_xleft=398.7,cl_ybottom=0.25,cl_xright=403,cl_ytop=0.30)
dev.off()

###########################################
###Apply Pareto Front functions for NHPP###
###########################################
## Focus on multiple testing units
source("func_Pareto_Front.R")
# Pareto Front optimization
dmv.nhpp.df_5_each_c=tp.eval.nhpp(ths=seq(0.1*365,1*365,by=5), 
                                  cs=seq(0,25,by=1),
                                  m1=0.0132,m0=0.009,
                                  pars.drs=waymo.bayes.we.altp.gibbs.obj$post.mat,
                                  mpd.test=0.2,mpd.use=0.2,
                                  dem.dur=2*365, nt=5)

tp.nhpp.pf.df_each_c=tp.dat(dmv.nhpp.df_5_each_c,
                            ths=seq(0.1*365,1*365,by=5),
                            cs=seq(0,25,by=1))

# Find the Pareto optimal test plans for a selected cut off based on consumer's risk
tp.nhpp.pf_nhpp=pf.cond.cr(tp.set=tp.nhpp.pf.df_each_c,cr.cut=0.13)

# Plot the Pareto optimal test plans for the NHPP model
fig.paper(filename="NHPP_Pareto_Front_each_c",type="pdf")
plot.pf.nhpp(tp.nhpp.pf_nhpp[[1]])
dev.off()