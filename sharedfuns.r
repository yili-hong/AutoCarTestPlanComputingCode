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