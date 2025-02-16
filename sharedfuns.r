########################################################################################
fig.paper=function(filename, type="eps",width=7,height=7)
{
  if(Sys.info()[['sysname']]=="Windows")
  {
    dir.create("./figures", showWarnings = F)
  }

  if(type=="ps" | type=="eps")
  {
   postscript(file=paste("./figures/",filename,".",type,sep=""),width=width,
   height=height,horizontal=F)
  }

  if(type=="eps.sq")
  {
   postscript(file=paste("./figures/",filename,".sq.eps",sep=""),width=7,
   height=7,horizontal=F)
  }

  if(type=="eps.rc")
  {
   postscript(file=paste("./figures/",filename,".rc.eps",sep=""),width=11,
   height=7,horizontal=T)
  }


  if(type=="cairo.ps.sq")
  {
    library(grDevices)
    cairo_ps(filename=paste("./figures/",filename,".sq.eps",sep=""),width=7,
   height=7)
  }

  if(type=="cairo.ps.rc")
  {
    library(grDevices)
    cairo_ps(filename=paste("./figures/",filename,".rc.eps",sep=""), width=11,
   height=7)
  }


  if(type=="pdf")
  {
   pdf(file=paste("./figures/",filename,".",type,sep=""),width=width,
   height=height)
  }
  
  if(type=="wmf")
  {
   win.metafile(filename = paste("./figures/",filename,".wmf",sep=""), 
   width = width, height = height, pointsize = 12, restoreConsole = TRUE)
  }

  if(type=="emf")
  {
   library(devEMF)
   emf(file = paste("./figures/",filename,".emf",sep=""), 
   width = width, height = height, pointsize = 12)
  }
  
  if(type=="tiff")
  {
   tiff(filename = paste("./figures/",filename,".tiff",sep=""),
   width = width, height = height, units = "in",
   pointsize = 12, res = 600, compression = "lzw")
  }

  if(type=="png")
  {
   png(filename = paste("./figures/",filename,".png",sep=""),
   width = width, height = height, units = "in",
   pointsize = 12, res = 600)
  }
  
  if(type=="jpeg")
  {
   jpeg(filename = paste("./figures/",filename,".jpeg",sep=""),
     width = width, height = height, units = "in", pointsize = 12,
     quality = 100,
     bg = "white", res = 360)

  }
}
########################################################################################
#####cdf of sev#####
psev=function(z)
{
  1-exp(-exp(z))
}
########################################################################################
#####density of sev#####
dsev=function(z)
{
  exp(z-exp(z))
}
########################################################################################
#####different parameterization of the Weibull model#####
weibull.fun.altp=function(tt, alpha, mu, sigma)
{
  zz=(log(tt)-mu)/sigma
  
  res=alpha*psev(zz)
  return(res)
}
########################################################################################
count.dat.setup=function(event.dat, mile.dat, int.tab)
{
  vins=unique(mile.dat[,"VIN"])
  mid=int.tab[,"MonthID"]
  
  Non.event.mat=NULL  
  Event.mat=NULL
  
  for(i in 1:length(vins))
  {
    vin=vins[i]
    yy=mile.dat[mile.dat[,"VIN"]==vin,]
    #browser()      
    ww=yy[mid]
    names(ww)=NULL
    yy1=data.frame(VIN=vin, Start=int.tab[,"start.day"], End=int.tab[,"end.day"], Mile=as.vector(t(ww)))
    Non.event.mat=rbind(Non.event.mat, yy1)
    
    if(vin%in%event.dat[,"VIN"])
    {
      xx=event.dat[event.dat[,"VIN"]==vin,]
      mm=t(yy[xx[,"MonthID"]])
      rownames(mm)=NULL
      xx1=data.frame(VIN=xx[,"VIN"], Start=xx[,"Days"]-0.001, End=xx[,"Days"], Mile=as.vector(mm))   
    }else{
      xx1=NULL
    }
    
    Event.mat=rbind(Event.mat, xx1)
    
  }
  
  Non.event.mat=Non.event.mat[Non.event.mat[,"Mile"]!=0,]
  row.names(Non.event.mat)=NULL
  
  Non.event.mat=data.frame(Non.event.mat, wts=1)
  Event.mat=data.frame(Event.mat, wts=1)
  
  res=list(Non.event.mat=Non.event.mat, Event.mat=Event.mat)
  
  return(res)
  
}
########################################################################################
data.setup=function(dat, manufacturer)
{
  event.dat=dat$event.dat
  mile.dat=dat$mile.dat
  initial.time=dat$initial.time
  time.int.tab=dat$time.int.tab
  
  event.dat=event.dat[event.dat[,"Manufacturer"]==manufacturer,]
  event.dat=event.dat[order(event.dat[,"VIN"], event.dat[,"Date"]),]
  event.dat=cbind(event.dat, Days=as.numeric(event.dat[,"Date"]-initial.time))
  
  time.int.tab=cbind(time.int.tab, start.day=as.numeric(time.int.tab[,"start.time"]-initial.time), end.day=as.numeric(time.int.tab[,"end.time"]-initial.time))
  
  mile.dat=mile.dat[mile.dat[,"Manufacturer"]==manufacturer,]
  monthIDs=time.int.tab[,"MonthID"]
  
  tmp=mile.dat[, monthIDs]/1000
  tmp1=sweep(tmp, 2, time.int.tab[,"len"],"/")
  
  mile.dat[, monthIDs]=tmp1
  mile.dat=mile.dat[rowSums(tmp)>0,]
  row.names(mile.dat)=NULL
  
  print("All event VINs contained in mile VINs?")
  print(all(unique(event.dat[,"VIN"]) %in% mile.dat[,"VIN"]))
  
  #browser()
  #aa=event.dat[!(event.dat[,"VIN"] %in% mile.dat[,"VIN"]),]
  #table(aa[,"VIN"])
  
  event.dat=event.dat[event.dat[,"VIN"] %in% mile.dat[,"VIN"],]
  rownames(event.dat)=NULL
  print(all(unique(event.dat[,"VIN"]) %in% mile.dat[,"VIN"]))
  
  
  count.obj=count.dat.setup(event.dat=event.dat, mile.dat=mile.dat, int.tab=time.int.tab)
  
  #browser()
  
  res=list(event.dat=event.dat, mile.dat=mile.dat, initial.time=initial.time, time.int.tab=time.int.tab, count.obj=count.obj)
  
  return(res)
  
}
########################################################################################
minus.loglik.weibull.altp=function(dat, pars)
{
  
  alpha=exp(pars[1])
  mu=pars[2] 
  sigma=exp(pars[3])
  
  #browser()
  
  Non.event.mat=dat$Non.event.mat
  Event.mat=dat$Event.mat
  
  aa1=weibull.der.fun.altp(tt=Event.mat[,"End"], alpha=alpha, mu=mu, sigma=sigma)
  pp1=(-1)*sum(log(aa1*Event.mat[,"Mile"])*Event.mat[,"wts"])
  
  bb1=weibull.fun.altp(tt=Non.event.mat[,"Start"], alpha=alpha, mu=mu, sigma=sigma)
  bb2=weibull.fun.altp(tt=Non.event.mat[,"End"], alpha=alpha, mu=mu, sigma=sigma)
  
  #browser()
  
  dd=bb2-bb1
  dd[dd==0]=max(dd)
  
  pp2=sum(Non.event.mat[,"Mile"]*dd*Non.event.mat[,"wts"])
  
  res=pp1+pp2
  
  return(res)
  
}
########################################################################################
lifetime.mle=function(dat, minusloglik, starts, method = "BFGS",hessian = TRUE,...)
{
  call=match.call()
  f = function(p) {
    minusloglik(dat,p) 
  }
  oout = optim(starts, f, method = method, hessian = hessian,...)#,control=list(trace=T))
  coef = oout$par
  #browser()
  if(hessian)
  {
    vcov =solve(oout$hessian)
  }else{
    vcov=NULL
  }
  min = oout$value
  invisible(list(call = call, coef = coef,vcov = vcov, min = min,dat=dat,minusloglik=minusloglik))
}
########################################################################################
Bayes.wei.altp.fit=function(dat, init, B=11000, B.sample=1000, b.sample=1, B.adapt = 1000, leap.eps=1, MCMC.method="Gibbs", prior.inf=NULL)
{
  
  cpp_obj=recurr.auto.car.dat.to.cpp(obj=dat, prior.inf=prior.inf)
  
  fun.run.start=proc.time()
  
  nuts_ctrl_list=list(int_M_adapt = B.adapt, num_delta = 0.5, int_max_treedepth = 10, num_eps = leap.eps)  
  
  if(MCMC.method=="NUTS")
  {
    
    #browser()
    
    cpp_obj[["vec_M_diag"]]=rep(1, length(init))
    
    post.mat=Nuts_wei_altp(vec_theta=init, int_iter=B+B.adapt, dat_list=cpp_obj, nuts_ctrl_list=nuts_ctrl_list)
    post.mat=post.mat[(B.adapt+1):(B.adapt+B),]
  }
  
  gibbs.obj=NULL
  
  if(MCMC.method=="Gibbs")
  {
    #browser()
    
    
    gibbs.obj=Gibbs_wei_altp(vec_theta=init, dat_list=cpp_obj, vec_scale=rep(1, length(init)), int_iter=B) 
    
    post.mat=gibbs.obj$vmat
    
    post.mat=post.mat[floor(seq(B.sample,B, by=b.sample)), ]
    
  }
  
  fun.run.end=proc.time()
  
  cat("MCMC running time: ", as.vector(fun.run.end-fun.run.start)[1:3], "\n")
  
  #######
  plot(post.mat[,ncol(post.mat)], type="l", xlab="Iteration", ylab="MCMC Chain", main="MCMC Chain Quick Check")
  
  res=list(obj=dat, cpp.obj=cpp_obj, post.mat=post.mat, gibbs.obj=gibbs.obj, nuts_ctrl_list=nuts_ctrl_list)   
  
  return(res)
  
}
########################################################################################
weibull.der.fun.altp=function(tt, alpha, mu, sigma)
{
  zz=(log(tt)-mu)/sigma
  tmp=dsev(zz)/(sigma*tt)   
  res=alpha*tmp
  
  return(res)
}
########################################################################################
recurr.auto.car.dat.to.cpp=function(obj, prior.inf=NULL)
{
  Non.event.mat=obj$Non.event.mat
  Event.mat=obj$Event.mat
  
  n1=nrow(Non.event.mat)
  n2=nrow(Event.mat)
  
  vin=c(Non.event.mat[,"VIN"], Event.mat[,"VIN"])
  tmp=as.numeric(as.factor(vin))
  
  Non.event.mat[,"VIN"]=tmp[1:n1]
  Event.mat[,"VIN"]=tmp[(n1+1):(n1+n2)]
  
  Non_event_mat=as.matrix(Non.event.mat)
  colnames(Non_event_mat)=NULL
  
  Event_mat=as.matrix(Event.mat)
  colnames(Event_mat)=NULL    
  
  #fix zeros.
  Non_event_mat[Non_event_mat[,2]==0,2]=1e-5 
  
  res=list(Non_event_mat=Non_event_mat, Event_mat=Event_mat, prior_inf=prior.inf)
  
  return(res)
  
}
########################################################################################