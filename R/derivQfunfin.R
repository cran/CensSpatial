derivQfun=function(est,fix.nugget=TRUE){

  if(class(est)!="SAEMSpatialCens") stop("an object of the class SAEMSpatialCens must be provided")
  if(!is.logical(fix.nugget)) stop("fix.nugget must be TRUE or FALSE")

  out=derivQ(est=est,fix.nugget=fix.nugget)

  return(list(Qlogvalue=out$Qlogvalue,gradQ=out$gradQ,HQ=out$HQ,QI=out$QI,Sigma=out$Sigma))

}
