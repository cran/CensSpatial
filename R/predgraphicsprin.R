
#########################################################################
############8. Prediction graphics for SAEM Algortihm##########################
##########################################################################

predgraph = function(xpred=NULL,grid1,est,points=TRUE,obspoints=1:sum(est$cc==0),colors=terrain.colors(100),sdgraph=TRUE,xlab="X Coord",ylab="Y Coord",main1="Predicted response", main2="Standard deviation predicted",xlim,ylim){

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
      xpred=xpred
    }



  }



  predgrap=predictionsaem(xpred=xpred,coordspred=coordspred,est=est)
  pred=predgrap

grid1$z1=pred$prediction

t=levelplot(z1~x*y,grid1,cuts = 30,col.regions=colors,main=main1,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
print(t)
if(points==TRUE){
trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(pred$coordsobs[est$cc==0,][obspoints,],pch=19,col=1,cex=0.5)
ltext(pred$coordsobs[est$cc==0,][obspoints,], labels =round(est$uy[est$cc==0],3), cex=0.5,pos=3)
trellis.unfocus()
}

on.exit(t)

  if(sdgraph==TRUE){
    grid1$z2=pred$sdpred
    dev21=dev.new(noRStudioGD=TRUE)
    on.exit(dev21)
sdp=levelplot(z2~x*y,grid1,cuts = 30,col.regions=colors,main=main2,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
print(sdp)
on.exit(sdp)
  }
  if(sdgraph==FALSE){
    a1=data.frame(predgrap$prediction,predgrap$coordspred)
    return(invisible(list(datapred=a1)))
  }

  if(sdgraph==TRUE){
    a1=data.frame(predgrap$prediction,predgrap$coordspred)
    a2=data.frame(predgrap$sdpred,predgrap$coordspred)
    return(invisible(list(datapred=a1,datasdpred=a2)))
  }

}

