
#########################################################################
############8. Prediction graphics for SAEM Algortihm##########################
##########################################################################

predgraph = function(xpred=NULL,grid1,est,points=T,obspoints=1:length(est$cc==0),colors=terrain.colors(100),sdgraph=T,legend.args.pred=list(text='Predicted response', side=4, font=2, line=2.5, cex=0.8),
                   legend.args.sdpred=list(text='Standard deviation predicted', side=4, font=2, line=2.5, cex=0.8)){

  trend=est$trend
  
  coords=est$coords;uy=est$uy;cc=est$cc;tau2=est$tau2;kappa=est$kappa;type=est$type

 
    colnames(coords)=c("x","y")
  coordspred=as.matrix(grid1)
  nobsest=length(uy)
  npred1=length(grid1[,1])
  nobs=length(cc==0)
  ncens=length(cc==1)
  interaction=coordspred[,1]*coordspred[,2]

  if(trend=="cte"){
    xpred=as.matrix(rep(1,npred1))
  }

  if(trend=="1st"){
    xpred=as.matrix(cbind(1,coordspred))
  }

  if(trend=="2nd"){
    xpred=as.matrix(cbind(1,coordspred,(coordspred)^2,interaction))
  }

  if(trend=="other"){
    if(is.null(xpred)){
      stop("object of the class SAEMSpatialCens was calculated with 
           trend= other, specify the xpred matrix")
    }
    else{
      xpred==xpred
    }



  }



  predgrap=predictionsaem(xpred=xpred,coordspred=coordspred,est=est)
  pred=predgrap
  ymin=xmin=min(pred$coordsobs[,1])
  ymax=xmax=max(pred$coordsobs[,1])
  r=raster(ncol=50,nrow=50,vals=predgrap$prediction)
  extent(r) <- c(xmin,xmax,ymin,ymax)
  r.range <- c(round(minValue(r),digits=2), round(maxValue(r),digits=2))




  plot(r,col=colors,asp=1,axes=T,xlab= "X coord",
       ylab="Y coord",axis.args=list(at=round(seq(r.range[1],
  r.range[2], length=10),digits=2),
  labels=round(seq(r.range[1], r.range[2], length=10),digits=2),
   cex.axis=0.6),legend.args=legend.args.pred)

  if(points==T){
    points(pred$coordsobs[est$cc==0,][obspoints,],pch=20)
    text(pred$coordsobs[est$cc==0,][obspoints,], labels = round(est$uy[est$cc==0][obspoints],digits=2), cex=0.5,pos=3)
  }

  if(sdgraph==T){
    ymin=xmin=min(pred$coordsobs[,1])
    ymax=xmax=max(pred$coordsobs[,1])
    r=raster(ncol=50,nrow=50,vals=predgrap$sdpred)
    extent(r) <- c(xmin,xmax,ymin,ymax)
    r.range <- c(round(minValue(r),digits=2), round(maxValue(r),digits=2))
    dev.new(noRStudioGD=TRUE)
    plot(r,col=colors,asp=1,axes=T,xlab= "X coord", ylab="Y coord",axis.args=list(at=round(seq(r.range[1], r.range[2], length=10),digits=2),
                                                                                  labels=round(seq(r.range[1], r.range[2], length=10),digits=2),
                                                                                  cex.axis=0.6),
         legend.args=legend.args.sdpred)
  }
  if(sdgraph==F){
    a1=data.frame(predgrap$prediction,predgrap$coordspred)
    return(invisible(list(datapred=a1)))
  }

  if(sdgraph==T){
    a1=data.frame(predgrap$prediction,predgrap$coordspred)
    a2=data.frame(predgrap$sdpred,predgrap$coordspred)
    return(invisible(list(datapred=a1,datasdpred=a2)))
  }

}
