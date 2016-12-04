predgraphics=function(xpred=NULL,grid1,est,points=T, obspoints=1:sum(est$cc==0),colors=terrain.colors(100),sdgraph=T,xlab="X Coord",ylab="Y Coord",
main1="Predicted response", main2="Standard deviation predicted",xlim=c(min(est$coords[,1]),max(est$coords[,1])),ylim=c(min(est$coords[,2]),max(est$coords[,2]))){

  if(class(est)!="SAEMSpatialCens") stop("an object of the class SAEMSpatialCens must be provided")
  #if (trend!='cte' & trend!='1st' & trend!="2nd" & trend!="other") stop('trend must be cte, 1st, or 2nd')

  if(!is.null(xpred)){
    if(!is.numeric(xpred) & !is.data.frame(xpred)) stop("xpred must be a numeric matrix or data.frame")
    if(!is.matrix(xpred)) xpred=as.matrix(xpred)
}

  if(!is.data.frame(grid1)) stop("grid1 must be a data.frame")
  if(nrow(grid1)!=nrow(xpred)) stop("all values of xpred for the specified coordinates must be specified")
  if(!is.character(xlab)) stop("Invalid value for xlab, must be a character")
  if(!is.character(ylab)) stop("Invalid value for ylab, must be a character")
  if(!is.character(main1)) stop("Invalid value for main1, must be a character")
  if(!is.character(main2)) stop("Invalid value for main2, must be a character")
  if(length(xlim)!=2) stop("Invalid value for xlim")
  if(length(ylim)!=2) stop("Invalid value for ylim")

if(!is.logical(points)) stop("points must be TRUE or FALSE")
  if(length(which(obspoints<0))>0) stop("Correct observations must be specified for the prediction graph")
  div=obspoints%%1
  if(length(which(div>0))>0) stop("Correct observations must be specified for the prediction graph")
  if(!is.logical(sdgraph)) stop("sdgraph must be TRUE or FALSE")

  out=predgraph(xpred=xpred,grid1=grid1,est=est,points=points,obspoints=obspoints,sdgraph=sdgraph,xlab=xlab,ylab=ylab,main1=main1,main2=main2,xlim=xlim,ylim=ylim,colors=colors)

  return(out)

}
