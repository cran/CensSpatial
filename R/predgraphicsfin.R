predgraphics=function(xpred=NULL,grid1,est,points=T, obspoints=1:length(est$cc==0),colors=terrain.colors(100),sdgraph=T,legend.args.pred=list(text='Predicted response', side=4, font=2, line=2.5, cex=0.8),
                      legend.args.sdpred=list(text='Standard deviation predicted', side=4, font=2, line=2.5, cex=0.8)){

  if(class()!="SAEMSpatialCens") stop("an object of the class SAEMSpatialCens must be provided")
  #if (trend!='cte' & trend!='1st' & trend!="2nd" & trend!="other") stop('trend must be cte, 1st, or 2nd')

  if(!is.null(xpred)){

    if(!is.numeric(xpred) & !is.data.frame(xpred)) stop("xpred must be a numeric matrix or data.frame")
    if(!is.numeric(xpred)) xpred=as.matrix(xpred)

    if(!is.numeric(grid1) & !is.data.frame(grid1)) stop("xpred must be a numeric matrix or data.frame")
    if(!is.numeric(grid1)) grid1=as.matrix(grid1)
    if(nrow(grid1)!=nrow(xpred)) stop("all values of xpred for the specified coordinates must be specified")
  }

  if(!is.logical(points)) stop("points must be TRUE or FALSE")
  if(length(which(obspoints<0))>0) stop("Correct observations must be specified for the prediction graph")
  div=obspoints%%1
  if(length(which(div>0))>0) stop("Correct observations must be specified for the prediction graph")
  if(!is.logical(sdgraph)) stop("sdgraph must be TRUE or FALSE")

  out=predgraph(xpred=xpred,grid1=grid1,est=est,points=points,obspoints=obspoints,sdgraph=sdgraph)

  return(out)

}
