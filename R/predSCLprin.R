###############7.Prediction for the SAEM algorithm###################

predictionsaem=function(xpred,coordspred,est){
  xobs=as.matrix(est$X)
  xpred=as.matrix(xpred)
  colnames(est$coords)=c("x","y")
  colnames(coordspred)=c("x","y")
  coords=rbind(est$coords,coordspred)
  n=dim(xobs)[1]+ dim(xpred)[1]
  npred=dim(xpred)[1]
  s=dim(xobs)[2]
  nobs=dim(xobs)[1]
  pred.row= (nobs+1):(nobs+npred)
  beta=est$theta[1:s]
  phi=est$theta[(s+2)]
  sigma2=est$theta[(s+1)]
  tau2=est$nugget
  kappa=est$kappa
  type=est$type
  z=est$uy

  Sigma=varcov.spatial(coords,cov.model=type,cov.pars=c(sigma2,phi),nugget=tau2,kappa=kappa)$varcov
  pred=rep(0,n)
  pred[pred.row]=1
  prediction=(xpred%*%beta)+ ( (Sigma[pred==1,pred==0]%*%solve(Sigma[pred==0,pred==0]))%*%(z-(xobs%*%beta)))
  S=Sigma[pred==1,pred==1]-(Sigma[pred==1,pred==0]%*%solve(Sigma[pred==0,pred==0])%*%Sigma[pred==0,pred==1])
  sdpred=sqrt(diag(S))
  coordspred=coords[pred==1,]
  coordsobs=coords[pred==0,]
  return(list(prediction=prediction,indpred=pred,sdpred=sdpred,coordspred=coordspred,coordsobs=coordsobs))
}
