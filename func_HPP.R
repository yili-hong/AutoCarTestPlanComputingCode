########################################################################################
crit.hpp.intensity=function(test.dur,c,m1,m0,pars.drs,mpd.test,mpd.use,ndm=6)
{ 
  nmc=dim(pars.drs)[1] #1001
  hist.time<- 2*365 #we freeze lamda at 2 years
  lambda.hpp=rep(0,nmc)
  for(i in 1:nmc)
  {
    lambda.hpp[i]=exp(pars.drs[i,1])*(dsev((log(hist.time)-pars.drs[i,2])/exp(pars.drs[i,3]))/(exp(pars.drs[i,3])*hist.time))
  }
  # prob. of the test is passed: no more than c failures have been observed
  f1=ppois(c,mpd.test*lambda.hpp*test.dur)   #where m=mpd.test*lambda.hpp under the HPP model
  ## prob. of the test will NOT be passed
  f2=1-f1		# prob. of the test is failure
  pt=sum(f1)
  ## prob. of having a successful test (acceptance rate)
  # MCMC
  accept.prob=round(pt/nmc,ndm) # the estimator we purse
  if(sum(f2)==0)
  {
    prod.risk=0
  }else
  {
    prod.risk=round(sum(f2[(mpd.use*lambda.hpp)<m0])/sum(f2),ndm)  # producer's risk
  }
  if(sum(mpd.use*lambda.hpp>m1)==0)
  {
    cons.risk=0
  }else
  {
    cons.risk=round(sum(f1[(mpd.use*lambda.hpp)>m1])/pt,ndm)	# consumer's risk
  }
  return(c(cons.risk,prod.risk,accept.prob,test.dur))
}

########################################################################################
## User-defined function to evaluate test plans
## Explore a set of test plans over a range of test hours and maximum allowable failures
########################################################################################
tp.eval.hpp=function(ths,cs,m1,m0,pars.drs,mpd.test,mpd.use)
{
  ncs=length(cs) 
  nths=length(ths)
  crit.vals=array(0,c(ncs,nths,4)) 
  for(i in 1:ncs)
  {
    for(j in 1:nths)
    {
      crit.vals[i,j,]=crit.hpp.intensity(ths[j],cs[i],m1,m0,pars.drs,mpd.test,mpd.use)
    }
  }
  return(crit.vals)
}

########################################################################################
### Visualization for the test plan###
# pre-set up for the change the default setting of color.legend 
########################################################################################
color.legend<- function(xl, yb, xr, yt, legend, rect.col, cex = 1, align = "lt", 
                        gradient = "x", ...)
{
  oldcex <- par("cex")
  par(xpd = TRUE, cex = cex)
  gradient.rect(xl, yb, xr, yt, col = rect.col, nslices = length(rect.col), 
                gradient = gradient, border = NA)
  if (gradient == "x") {
    xsqueeze <- (xr - xl)/(2 * length(rect.col))
    textx <- seq(xl + xsqueeze, xr - xsqueeze, length.out = length(legend))
    if (match(align, "rb", 0)) {
      texty <- yb - 0.2 * strheight("O")
      textadj <- c(0.5, 1)
    }
    else {
      texty <- yt + 0.2 * strheight("O")
      textadj <- c(0.5, 0)
    }
  }
  else {
    ysqueeze <- (yt - yb)/(2 * length(rect.col))
    texty <- seq(yb + ysqueeze, yt - ysqueeze, length.out = length(legend))
    if (match(align, "rb", 0)) {
      textx <- xr + 0.2 * strwidth("O")
      textadj <- c(0, 0.5)
    }
    else {
      textx <- xl - 0.2 * strwidth("O")
      textadj <- c(1, 0.5)
    }
  }
  text(textx, texty, labels = legend, adj = textadj, ...)
  par(xpd = FALSE, cex = oldcex)
}

########################################################################################
### Test Plan for HPP###
########################################################################################
## (1) CR vs PR
plot.CR.PR=function(critvals,ths,cs,x_start=0.7,x_end=0.95,y_start=0,y_end=0.7,
                    t=c(500,1000,2000,3000,4000),pt=c(0,15,19,17,18),
                    legend_labels=c(expression(paste(tau, "=500")),expression(paste(tau, "=1000")),
                                    expression(paste(tau, "=2000")),expression(paste(tau, "=3000")),
                                    expression(paste(tau, "=4000"))),
                    col_labels= c("c=0","c=10","c=20","c=30","c=40","c=50"),
                    cl_xleft=0.704,cl_ybottom=0.18,cl_xright=0.708,cl_ytop=0.35){
  ncs=length(cs)
  nths=length(ths)
  cols=gray.colors(ncs, start = 0.1, end = 0.85)
  plot(c(x_start,x_end),c(y_start,y_end),
       xlab="Consumer's Risk",ylab="Producer's Risk",main="",type="n")
  for(i in 1:ncs)
  {
    lines(critvals[i,,1],critvals[i,,2],lty=1,col=cols[i])
    for (j in 1:length(t)){
      criteria_ths<- which(ths==t[j]) #set up typical values of total testing time for further evaluate
      points(critvals[i,criteria_ths,1],critvals[i,criteria_ths ,2],pch = pt[j])
    }
  }
  legend('bottomleft',
         legend_labels,
         pch = c(0, 15, 19, 17, 18),
         bty = 'n',
         cex = 0.80,
         y.intersp=0.8,x.intersp=0.8,text.width=0.5)
  
  color.legend(xl=cl_xleft,yb=cl_ybottom,xr=cl_xright,yt=cl_ytop,
              legend=col_labels,
              align="rb",
              rect.col=matrix(cols,nrow=length(cs),ncol=1),
              gradient = "y",
              cex=0.7)
}

## (2) TH vs AP
plot.TH.AP=function(critvals,ths,cs,x_adjust=500,cr=c(0.7,0.75,0.8,0.85,0.9),pt1=c(1:5),
                    pr=c(0.2,0.3,0.4,0.5,0.6),pt2=c(17,18,20,15,0),
                    l_x=3680,l_y=1.05,
                    legend_labels=c('CR=0.05','CR=0.1','CR=0.15','CR=0.2','CR=0.25','PR=0.05','PR=0.1','PR=0.15','PR=0.2','PR=0.25'),
                    width=570,col_labels= c("c=0","c=10","c=20","c=30","c=40","c=50"),
                    cl_xleft=5100.97,cl_ybottom=0.30,cl_xright=5180,cl_ytop=0.60){
  ncs=length(cs)
  nths=length(ths)
  cols=gray.colors(ncs, start = 0.1, end = 0.85)
  plot(critvals[1,,4],critvals[1,,3],xlim=c(min(ths),max(ths)+x_adjust),ylim=c(0,1),xlab="Total Vehicle Days for Testing",ylab="Acceptance Probability",main="",type="n")
  for(i in 1:ncs)
  {
    lines(critvals[i,,4],critvals[i,,3],lty=1,col=cols[i])
    for (j in 1:length(cr)){
      # given target CR value, find corresponding values for total vehicle days
      ths_target_cr<- approx(x = critvals[i,,1], y = critvals[i,,4], xout = cr[j])$y
      # returns back to the AP vs. total testing time
      pr_target<- approx(x = critvals[i,,4], y = critvals[i,,3], xout = ths_target_cr)$y 
      points(ths_target_cr,pr_target,pch = pt1[j]) ## add points to the above lines
    }
    for (k in 1:length(pr)){
      # given target PR value, find corresponding values for total vehicle days
      ths_target_pr<- approx(x = critvals[i,,2], y = critvals[i,,4], xout = pr[k])$y
      # returns back to the AP vs. total testing time
      ap_target<- approx(x = critvals[i,,4], y = critvals[i,,3], xout = ths_target_pr)$y 
      points(ths_target_pr,ap_target,pch = pt2[k]) ## add points to the above lines
    }
  }
  legend(x=l_x,y=l_y,
         legend=legend_labels,
         pch = c(1,2,3,4,5,18, 17, 20, 15, 0),
         bty = 'n',
         cex = 0.8,
         inset = c(0, -0.02),
         ncol=2,
         y.intersp=0.8,
         x.intersp=0.8,
         text.width=width) 
 
  color.legend(xl=cl_xleft,yb=cl_ybottom,xr=cl_xright,yt=cl_ytop,
               legend=col_labels,
               align="rb",
               rect.col=matrix(cols,nrow=length(cs),ncol=1),
               gradient = "y",
               cex=0.7)
}

## (3) TH vs CR (Bottom Left)
plot.TH.CR=function(critvals,ths,cs,ystart=0,yend=0.5,
                    pr=c(0.05,0.10,0.15,0.2,0.25),pt=c(18,17,20,15,0),
                    l_x=2200,l_y=0.52,
                    legend_labels=c('PR=0.05','PR=0.1','PR=0.15','PR=0.2','PR=0.25'),
                    col_labels= c("c=0","c=10","c=20","c=30","c=40","c=50"),
                    cl_xleft=2300.97,cl_ybottom=0.23,cl_xright=2350,cl_ytop=0.38){
  ncs=length(cs)
  nths=length(ths)
  cols=gray.colors(ncs, start = 0.1, end = 0.85)
  plot(critvals[1,,4],critvals[1,,1],
       xlim=c(min(ths),max(ths)),
       ylim=c(ystart,yend),
       xlab="Total Vehicle Days for Testing",ylab="Consumer's Risk",main="",type="n")
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
         legend_labels,
         pch = c(18, 17, 20, 15, 0),
         bty = 'n',
         cex = 0.8,
         y.intersp=0.8,
         x.intersp=0.8,
         text.width=40,
         inset=c(0.05,0))

  color.legend(xl=cl_xleft,yb=cl_ybottom,xr=cl_xright,yt=cl_ytop, 
               legend=col_labels,
               align="rb",
               rect.col=matrix(cols,nrow=length(cs),ncol=1),
               gradient = "y",
               cex=0.7)
}

## (4) TH vs PR (Bottom Right)
plot.TH.PR=function(critvals,ths,cs,ystart=0,yend=0.3,
                    cr=c(0.7,0.75,0.8,0.85,0.9),pt=c(1:5),
                    x_adjust,l_x=2170,l_y=0.52, 
                    legend_labels=c('CR=0.07','CR=0.09','CR=0.11','CR=0.13','CR=0.15'),
                    cl_xleft=2300.97,cl_ybottom=0.07,cl_xright=2350,cl_ytop=0.17,
                    col_labels= c("c=0","c=10","c=20","c=30","c=40","c=50")){
  ncs=length(cs)
  nths=length(ths)
  cols=gray.colors(ncs, start = 0.1, end = 0.85)
  plot(critvals[1,,4],
       critvals[1,,2],
       xlim=c(min(ths),max(ths)+x_adjust),
       ylim=c(ystart,yend),
       xlab="Total Vehicle Days for Testing",ylab="Producer's Risk",main="",type="n")
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
         legend_labels,
         pch = c(1, 2, 3, 4, 5),
         bty = 'n',
         cex = 0.8,
         y.intersp=0.8,
         x.intersp=0.8,
         text.width=30,
         inset=c(0.05,0))

  color.legend(xl=cl_xleft,yb=cl_ybottom,xr=cl_xright,yt=cl_ytop, 
               legend=col_labels,
               align="rb",
               rect.col=matrix(cols,nrow=length(cs),ncol=1),
               gradient = "y",
               cex=0.7)
}
########################################################################################