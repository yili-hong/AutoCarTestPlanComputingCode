########################################################################################
##### Function to calculate the criteria based on #####
# NHPP with Weibull model                      
# using average failure intensity as reliability demonstration metric  
########################################################################################
crit.nhpp.weibull.rel=function(test.dur,c,m1,m0,pars.drs,mpd.test,mpd.use,dem.dur,nt,ndm=6)
{
  nmc=dim(pars.drs)[1] 
  hist.time<- 2*365 #2 years
  
  Lambda_test_st=Lambda_test_end=a_Lambda_test=Lambda_dem_st=Lambda_dem_end=a_Lambda_dem=rep(0,nmc)
  f1=f2=f1s=f2s=rep(0,length(nmc)) #only one point
  for(i in 1:nmc) 
  { 
    # CBIF for the testing period
    Lambda_test_st[i] = psev((log(hist.time)-pars.drs[i,2])/exp(pars.drs[i,3]))*exp(pars.drs[i,1])
    Lambda_test_end[i] = psev((log(hist.time+test.dur)-pars.drs[i,2])/exp(pars.drs[i,3]))*exp(pars.drs[i,1])
    a_Lambda_test[i] = (Lambda_test_end[i]-Lambda_test_st[i])/test.dur 
   
    ppois(c,mpd.test*a_Lambda_test[i]*test.dur*nt)
    #prob. of the test is passed if y<= c
    f1[i]=ppois(c,mpd.test*a_Lambda_test[i]*test.dur*nt)
    #prob. of failure for each row if y>c
    f2[i]=1-f1[i]		
    
    # CBIF for the demonstration period
    Lambda_dem_st[i] = psev((log(hist.time)-pars.drs[i,2])/exp(pars.drs[i,3]))*exp(pars.drs[i,1])
    Lambda_dem_end[i] = psev((log(hist.time+dem.dur)-pars.drs[i,2])/exp(pars.drs[i,3]))*exp(pars.drs[i,1])
    a_Lambda_dem[i] = (Lambda_dem_end[i]-Lambda_dem_st[i])/dem.dur 
    # Note that this is the /mileage-adjusted cumulative intensity/ for the form of Weibull model
    
    # f1s: joint prob. of the test is passed and the prob. of the posterior consumer's risk
    # f2s: joint prob. of the test of failed and the prob. of the posterior producer's risk
    # as.numeric: posterior consumer's risk
    f1s[i]=f1[i]*as.numeric((mpd.use*a_Lambda_dem[i])>m1)  
    # as.numeric: posterior producer's risk
    f2s[i]=f2[i]*as.numeric((mpd.use*a_Lambda_dem[i])<m0)
  }
  
  pt=sum(f1)
  accept.prob=round(pt/nmc,ndm) #averaged pass probability, round to 6 decimal digit
  if(sum(f2)==0)
  {
    prod.risk=0
  }else
  {
    prod.risk=round(sum(f2s)/sum(f2),ndm)  # producer's risk-> f2s:joint prob.; f2: prob of failed
  }
  if(sum((mpd.use*a_Lambda_dem[i])>m1)==0)
  {
    cons.risk=0
  }else
  {
    cons.risk=round(sum(f1s)/sum(f1),ndm)	# consumer's risk+ 
  }
  return(c(cons.risk,prod.risk,accept.prob,test.dur))
}

########################################################################################  
##### User-defined Function to evaluate test plans #####
# explore a set of test plans over a range of test hours and maximum allowable failures
########################################################################################
tp.eval.nhpp=function(ths,cs,m1,m0,pars.drs,mpd.test,mpd.use,dem.dur,nt)
{
  ncs=length(cs) 
  nths=length(ths) 
  crit.vals=array(0,c(ncs,nths,4)) #(t*,c) which represent the test plan, where t* is the testing period
  for(i in 1:ncs)
  {
    for(j in 1:nths)
    {
      crit.vals[i,j,]=crit.nhpp.weibull.rel(ths[j],cs[i],m1,m0,pars.drs,mpd.test,mpd.use,dem.dur,nt)
    }
  }
  return(crit.vals)
}

########################################################################################
##### Test plan for NHPP #####
########################################################################################
# pre-set up for the change the default setting of color.legend 
source("color_legend.R")
## (1) CR vs PR
plot.CR.PR=function(critvals,ths,cs,x_start=0,x_end=0.5,y_start=0,y_end=0.8,
                    t=c(51.5,101.5,201.5,301.5),pt=c(0,15,19,17,18),
                    legend_labels=c(expression(paste(tau[t], "=50")),expression(paste(tau[t], "=100")),
                                    expression(paste(tau[t], "=200")),expression(paste(tau[t], "=300"))),
                    color_labels=c("c=0","c=2","c=4","c=6","c=8","c=10"),
                    cl_xleft=0.004,cl_ybottom=0.18,cl_xright=0.008,cl_ytop=0.35)
{
  ncs=length(cs)
  nths=length(ths)
  
  cols=gray.colors(ncs, start = 0.1, end = 0.85)
  plot(c(x_start,x_end),c(y_start,y_end),
       xlab="Consumer's Risk",ylab="Producer's Risk",main="",type="n")
  for(i in 1:ncs)
  {
    lines(critvals[i,,1],critvals[i,,2],lty=1,col=cols[i])
    for (j in 1:length(t)){
      # given target total vehicle days, find corresponding values for CR value
      ths_target_cr<- approx(x = critvals[i,,4], y = critvals[i,,1], xout = t[j])$y
      # returns back to the CR (given) vs. PR
      pr_target<- approx(x = critvals[i,,1], y = critvals[i,,2], xout = ths_target_cr)$y 
      points(ths_target_cr,pr_target,pch = pt[j])
    }
  }
  legend('bottomleft',
         legend_labels,
         pch = pt,
         bty = 'n',
         cex = 0.80,
         y.intersp=0.8,x.intersp=0.8,text.width=0.5)
  col_labels<- color_labels
  color.legend(xl=cl_xleft,yb=cl_ybottom,xr=cl_xright,yt=cl_ytop,
               legend=col_labels,
               align="rb",
               rect.col=matrix(cols,nrow=length(cs),ncol=1),
               gradient = "y",
               cex=0.7)
}

## (2) TH vs AP
plot.TH.AP=function(critvals,ths,cs,x_adjust=500,
                    cr=c(0.7,0.75,0.8,0.85,0.9),pt1=c(1:5),
                    pr=c(0.2,0.3,0.4,0.5,0.6),pt2=c(17,18,20,15,0),
                    l_x=360,l_y=1.0, 
                    legend_label=c('CR=0.05','CR=0.1','CR=0.15','CR=0.2','CR=0.25','PR=0.05','PR=0.1','PR=0.15','PR=0.2','PR=0.25'),
                    width=70, color_labels=c("c=0","c=10","c=20","c=30","c=40","c=50"),
                    cl_xleft=5100.97,cl_ybottom=0.30,cl_xright=5180,cl_ytop=0.60){
  ncs=length(cs)
  nths=length(ths)
  cols=gray.colors(ncs, start = 0.1, end = 0.85)
  plot(critvals[1,,4],critvals[1,,3],
       xlim=c(min(ths),max(ths)+x_adjust),
       ylim=c(0,1),xlab="Total Testing Period",
       ylab="Acceptance Probability",main="",type="n")
  for(i in 1:ncs)
  {
    lines(critvals[i,,4],critvals[i,,3],lty=1,col=cols[i])
    for (j in 1:5){
      # given target CR value, find corresponding values for total vehicle days
      ths_target_cr<- approx(x = critvals[i,,1], y = critvals[i,,4], xout = cr[j])$y
      # returns back to the AP vs. total testing time
      pr_target<- approx(x = critvals[i,,4], y = critvals[i,,3], xout = ths_target_cr)$y 
      points(ths_target_cr,pr_target,pch = pt1[j]) ##add points to the above lines
    }
    for (k in 1:5){
      # given target PR value, find corresponding values for total vehicle days
      ths_target_pr<- approx(x = critvals[i,,2], y = critvals[i,,4], xout = pr[k])$y
      # returns back to the AP vs. total testing time
      ap_target<- approx(x = critvals[i,,4], y = critvals[i,,3], xout = ths_target_pr)$y 
      points(ths_target_pr,ap_target,pch = pt2[k]) ## add points to the above lines
    }
  }
  legend(x=l_x,y=l_y,
         legend=legend_label,
         pch = c(1,2,3,4,5,18, 17, 20, 15, 0),
         bty = 'n',
         cex = 0.8,
         inset = c(0, -0.02),
         ncol=2,
         y.intersp=0.8,
         x.intersp=0.8,
         text.width=width) 
  col_labels<- color_labels
  color.legend(xl=cl_xleft,yb=cl_ybottom,xr=cl_xright,yt=cl_ytop,
               legend=col_labels,
               align="rb",
               rect.col=matrix(cols,nrow=length(cs),ncol=1),
               gradient = "y",
               cex=0.7)
}

# (3) TH vs CR (Bottom Left)
plot.TH.CR=function(critvals,ths,cs,ystart=0,yend=0.5,
                    pr=c(0.05,0.10,0.15,0.2,0.25),pt=c(18,17,20,15,0),
                    l_x=2200,l_y=0.52, legend_label=c('PR=0.05','PR=0.1','PR=0.15','PR=0.2','PR=0.25'),
                    color_labels=c("c=0","c=10","c=20","c=30","c=40","c=50"),
                    cl_xleft=2300.97,cl_ybottom=0.23,cl_xright=2350,cl_ytop=0.38){
  ncs=length(cs)
  nths=length(ths)
  cols=gray.colors(ncs, start = 0.1, end = 0.85)
  plot(critvals[1,,4],critvals[1,,1],
       xlim=c(min(ths),max(ths)),
       ylim=c(ystart,yend),xlab="Total Testing Period",ylab="Consumer's Risk",main="",type="n")
  for(i in 1:ncs)
  {
    lines(critvals[i,,4],critvals[i,,1],lty=1,col=cols[i])
    for (j in 1:length(pr)){
      # given target PR value, find corresponding values for total vehicle days
      ths_target<- approx(x = critvals[i,,2], y = critvals[i,,4], xout = pr[j])$y
      # returns back to the CR vs. total testing time
      cr_target<- approx(x = critvals[i,,4], y = critvals[i,,1], xout = ths_target)$y 
      points(ths_target,cr_target,pch = pt[j]) ## add points to the above lines
    }
  }
  legend(x=l_x,y=l_y,
         legend_label,
         pch = c(18, 17, 20, 15, 0),
         bty = 'n',
         cex = 0.8,
         y.intersp=0.8,
         x.intersp=0.8,
         text.width=40,
         inset=c(0.05,0))
  col_labels<- color_labels
  color.legend(xl=cl_xleft,yb=cl_ybottom,xr=cl_xright,yt=cl_ytop, 
               legend=col_labels,
               align="rb",
               rect.col=matrix(cols,nrow=length(cs),ncol=1),
               gradient = "y",
               cex=0.7)
}

## (4) TH vs PR (Bottom Right)
plot.TH.PR=function(critvals,ths,cs,x_adjust=100,ystart=0,yend=0.3,
                    cr=c(0.7,0.75,0.8,0.85,0.9),pt=c(1:5),
                    l_x=2170,l_y=0.52,legend_label=c('CR=0.05','CR=0.1','CR=0.15','CR=0.2','CR=0.25'),
                    color_labels=c("c=0","c=2","c=4","c=6","c=8","c=10"),
                    cl_xleft=2300.97,cl_ybottom=0.07,cl_xright=2350,cl_ytop=0.17){
  ncs=length(cs)
  nths=length(ths)
  cols=gray.colors(ncs, start = 0.1, end = 0.85)
  plot(critvals[1,,4],critvals[1,,2],
       xlim=c(min(ths),max(ths)+x_adjust),
       ylim=c(ystart,yend),
       xlab="Total Testing Period",
       ylab="Producer's Risk",main="",type="n")
  for(i in 1:ncs)
  {
    lines(critvals[i,,4],critvals[i,,2],lty=1,col=cols[i])
    for (j in 1:length(cr)){
      # given target CR value, find corresponding values for total vehicle days
      ths_target<- approx(x = critvals[i,,1], y = critvals[i,,4], xout = cr[j])$y
      # returns back to the PR vs. total testing time
      pr_target<- approx(x = critvals[i,,4], y = critvals[i,,2], xout = ths_target)$y 
      points(ths_target,pr_target,pch = pt[j]) ## add points to the above lines
    }
  }
  legend(x=l_x,y=l_y,
         legend_label,
         pch = c(1, 2, 3, 4, 5),
         bty = 'n',
         cex = 0.8,
         y.intersp=0.8,
         x.intersp=0.8,
         text.width=30,
         inset=c(0.05,0))
  col_labels<- color_labels
  color.legend(xl=cl_xleft,yb=cl_ybottom,xr=cl_xright,yt=cl_ytop, 
               legend=col_labels,
               align="rb",
               rect.col=matrix(cols,nrow=length(cs),ncol=1),
               gradient = "y",
               cex=0.7)
}
########################################################################################