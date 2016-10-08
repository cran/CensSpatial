### 9.Maximum Likelihood Expectation (Q function and its derivates)###




derivQ=function(est,fix.nugget=T){
  qlog=function(theta,coords,X,cov.model,kappa,uy,uyy){
    k=dim(X)[2]
    beta=theta[1:k]
    if(fix.nugget==T){
      tau2=est$nugget
    }
    else{
      tau2=theta[(k+3)]}
    sigma2=theta[(k+1)]
    phi=theta[(k+2)]
    Psi<-varcov.spatial(coords,cov.model=type,cov.pars=c(sigma2,phi),nugget=tau2,kappa=kappa)$varcov
    V1=solve(Psi)
    media=X%*%beta
    r=-0.5*(log(det(Psi))+tr(uyy%*%V1)-(2*uy%*%V1%*%media) + (t(media)%*%V1%*%media))
    return(r)
  }

  beta=est$beta
  sigma2=est$sigma2
  phi=est$phi
  tau2=est$nugget
  if(fix.nugget==T){
    theta=c(beta,sigma2,phi)
  }
  else{
    theta=c(beta,sigma2,phi,tau2)
  }

  type=est$type
  uy=est$uy
  uyy=est$uyy
  beta=est$beta
  X=as.matrix(est$X)
  coords=est$coords
  Psi<-varcov.spatial(coords,cov.model=type,cov.pars=c(sigma2,phi),nugget=tau2,kappa=kappa)$varcov
  r=qlog(theta,coords=coords,X=X,cov.model=type,kappa=kappa,uy=uy,uyy=uyy)
  s=grad(qlog,theta,coords=coords,X=X,cov.model="matern",kappa=kappa,uy=uy,uyy=uyy)
  Q=hessian(qlog,theta,coords=coords,X=X,cov.model="matern",kappa=kappa,uy=uy,uyy=uyy)
  Qinv=solve(-Q)

  return(list(Qlogvalue=r,gradQ=s,HQ=Q,QI=Qinv,Sigma=Psi))

}

