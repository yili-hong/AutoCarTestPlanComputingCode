# Computing Code for Test Planning of Autonomous Vehicles

This public GitHub repository contains the R Program for the paper "Planning Reliability Assurance Tests for Autonomous Vehicles Based on Disengagement Events Data" by Simin Zheng, Lu Lu, Yili Hong, and Jian Liu. The paper is published at IISE Transactions. 

This public GitHub repository contains 12 files, listed below.

(1) readme.md  

(2) data folder: contains CA DMV disengagement data for Waymo used in this paper. Specifically, there are three .csv files under the data folder:
     (i) event.dat.csv: contains the event times.
    (ii) mile.dat.csv: contains the mileage information. 
   (iii) time.int.tab.csv: contains time interval information. 

(3) script_data_MCMC.R: Main R scripts for generating the Waymo MCMC, saved as the output in equation (3) of the paper .

(4) waymo.bayes.we.altp.gibbs.obj: The posterior inputs used in equations (8) and (9) of the paper.

(5) sharedfuns.r: Contains functions for general-purpose use.

(6) functions.cpp: Rcpp functions used to generate the Waymo MCMC. Needs to be compiled before they can be used in R.

(7) func_HPP.R: Contains functions related to the HPP models used in this paper. 

(8) func_NHPP.R: Contains functions related to the NHPP models used in this paper.

(9) func_Pareto_Front.R: Contains functions for Pareto Front optimization.

(10) color_legend.R: Contains functions for generating visualizations. 

(11) script_HPP.R: Main R scripts for examples and Pareto Front optimal test plans for the HPP model.

(12) script_NHPP.R: Main R scripts for examples and Pareto Front optimal test plans for the NHPP model.

