predSCL=function(xpred,coordspred,est){
  if(!is.numeric(xpred) & !is.data.frame(xpred)) stop ("xpred must be a numeric matrix or data.frame")
  if(!is.numeric(coordspred) & !is.data.frame(coordspred)) stop ("coordspred must be a numeric matrix or data.frame")
  if (!is.matrix(xpred)) xpred=as.matrix(xpred)
  if (!is.matrix(coordspred)) as.matrix(coordspred)
  if(ncol(coordspred)!=2) stop("2D coordinates must be specified")
  if(nrow(xpred)!=nrow(coordspred)) stop("xpred does not have the same number of lines than coordspred")
  if(!inherits(est,'SAEMSpatialCens')) stop("an object of the class SAEMSpatialCens must be provided")

  if(sum(is.na(xpred)) > 0) stop("There are some NA values in xpred")
  if(sum(is.na(coordspred)) > 0) stop("There are some NA values in coordspred")


  out=predictionsaem(xpred=xpred,coordspred=coordspred,est=est)

  return(out)

}
