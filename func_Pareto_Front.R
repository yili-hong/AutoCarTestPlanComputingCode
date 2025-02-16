################################################################################
##### Pareto Front related functions#####
###User-function to identify the Pareto front###
# based on pr(producer's risk, to minimize), ap(acceptance probability, to maximize), and th (total test hour, to minimize)
checkon3c=function(newpt,curpf)
{
  ##each criteria for the new line is better than the previous line
  g1=newpt[1,2]<curpf[,2] #this is for pr (min)
  g2=newpt[1,3]>curpf[,3] #this is for ap (max)
  g3=newpt[1,4]<curpf[,4] #this is for total testing time (min)
  
  #includes equal signs for ge1,ge2, and ge3
  ge1=newpt[1,2]<=curpf[,2]
  ge2=newpt[1,3]>=curpf[,3]
  ge3=newpt[1,4]<=curpf[,4]
  
  ##each criteria for the new line is worse than the previous line
  l1=newpt[1,2]>curpf[,2]
  l2=newpt[1,3]<curpf[,3]
  l3=newpt[1,4]>curpf[,4]
  #includes equal signs
  le1=newpt[1,2]>=curpf[,2]
  le2=newpt[1,3]<=curpf[,3]
  le3=newpt[1,4]>=curpf[,4]
  
  ##the same 
  eq1=newpt[1,2]==curpf[,2]
  eq2=newpt[1,3]==curpf[,3]
  eq3=newpt[1,4]==curpf[,4]
  
  ##key parts are below
  cond1=(g1*ge2*ge3+g2*ge1*ge3+g3*ge1*ge2)==0 #TRUE or FALSE
  cond2=sum(l1*le2*le3+l2*le1*le3+l3*le1*le2+eq1*eq2*eq3) # cond2=1 means add the point
  
  newpf=matrix(curpf[cond1,],ncol=dim(curpf)[2]) #keep the original row
  if(cond2==0) #if the criteria is not worse than the previous one, then we will add this row
  {
    newpf=rbind(newpf,newpt)
  }
  return(newpf)
}

################################################################################

###Create a data frame for the criteria values for the explored range of test plans###
tp.dat=function(critvals,ths,cs)
{
  tp.df=cbind(rep(cs[1],length(ths)),critvals[1,,])
  for(i in 2:length(cs))
  {
    tp.df=rbind(tp.df,cbind(rep(cs[i],length(ths)),critvals[i,,]))
  }
  #tp.df=as.data.frame(tp.df)
  colnames(tp.df)=c("c","cons.risk","prod.risk","accept.prob","test.time")
  return(tp.df)
}

###Find the Pareto front based on producer's risk, acceptance probability and total test time
# for test plans satisfying the primary criterion based on the consumer's risk cut off###
pf.cond.cr=function(tp.set,cr.cut)
{
  #start time
  time_st<- proc.time()[3]
  tp.set2=tp.set[tp.set[,"cons.risk"]<cr.cut,]
  pfs=matrix(tp.set2[1,c(1,3,4,5,2)],nrow=1) 
  for(i in 2:dim(tp.set2)[1]) 
  {
    pfs=checkon3c(matrix(tp.set2[i,c(1,3,4,5,2)],nrow=1),pfs) 
    }
  pfs2=pfs[,c(1,5,2,3,4)]
  colnames(pfs2)=c("c","cons.risk","prod.risk","accept.prob","test.time")
  #end time
  time_end<- proc.time()[3]
  class(time_end)
  return(list(pfs2,time = time_end-time_st))
}

################################################################################

###Plot the Pareto Front Test Plans for the HPP model###
plot.pf=function(pf.dat)
{
  pf.set=as.data.frame(pf.dat)
  names(pf.set)=colnames(pf.dat)
  tt.range=max(pf.set$test.time)-min(pf.set$test.time) 
  if(tt.range==0)
  {
    pf.set$tt.sc=rep(1,dim(pf.set)[1])
  }else
  {
    pf.set$tt.sc=(pf.set$test.time-min(pf.set$test.time))/(tt.range) 
  }
  par(mai=c(1,1,0.5,1))
  plot(pf.set$c,pf.set$prod.risk,
       xlab="Maximum Allowed Number of Disengagement Events",
       ylab="Probability",
       ylim=c(0,1),main="",pch=16,cex=1.2,axes=F)
  axis(1,at=seq(0,50,by=5),labels=seq(0,50,by=5)) 
  axis(2,at=seq(0,1,by=0.2),labels=seq(0,1,by=0.2)) 
  axis(4,at=seq(0,1,by=0.2),labels=round(min(pf.set$test.time)+seq(0,1,by=0.2)*tt.range,0))
  mtext("Total Vehicle Days for Testing (Days)",side=4,line=3)
  points(pf.set$c,pf.set$accept.prob,pch=17,cex=1.2)
  points(pf.set$c,pf.set$tt.sc,pch=15,cex=1.2)
  legend("topleft",
         pch=c(16,17,15),
         bty="n",
         legend=c("Producer's Risk","Acceptance Probability","Total Vehicle Days for Testing"))
}

################################################################################

###Plot the Pareto Front Test Plans for the NHPP###
plot.pf.nhpp=function(pf.dat)
{
  pf.set=as.data.frame(pf.dat)
  names(pf.set)=colnames(pf.dat)
  tt.range=max(pf.set$test.time)-min(pf.set$test.time) 
  if(tt.range==0)
  {
    pf.set$tt.sc=rep(1,dim(pf.set)[1])
  }else
  {
    pf.set$tt.sc=(pf.set$test.time-min(pf.set$test.time))/(tt.range) 
  }
  par(mai=c(1,1,0.5,1))
  plot(pf.set$c,pf.set$prod.risk,
       xlab="Maximum Allowed Number of Disengagement Events",
       ylab="Probability",
       ylim=c(0,1),main="",pch=16,cex=1.2,axes=F)
  axis(1,at=seq(0,25,by=5),labels=seq(0,25,by=5)) 
  axis(2,at=seq(0,1,by=0.2),labels=seq(0,1,by=0.2)) 
  axis(4,at=seq(0,1,by=0.2),labels=round(min(pf.set$test.time)+seq(0,1,by=0.2)*tt.range,0))
  mtext("Total Test Time (Days)",side=4,line=3)
  points(pf.set$c,pf.set$accept.prob,pch=17,cex=1.2)
  points(pf.set$c,pf.set$tt.sc,pch=15,cex=1.2)
  legend("topleft",
         pch=c(16,17,15),
         bty="n",
         legend=c("Producer's Risk","Acceptance Probability","Total Testing Vehicle Days"))
}
################################################################################