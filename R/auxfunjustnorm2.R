
#####################################
######################################3
#######################################

################ALGORITMO SAEM (NORMAL CASE)############################################
################################################################################
## Funcao para maximizar usada no SAEM: Estima phi e tau!!!
################################################################################

FCi<-function(rhoG,media,sigmae,kappa,uu,yb,uyb,yyb,type.S,coords){
  phi<-rhoG[1]
  tau2<-rhoG[2]
  vero=numeric()

  y1 <- yb
  u1 <- uu
  uy1 <- uyb
  yy1 <- yyb

  #Cii<-MatSpatial(MCoord,sigmae, phi, tau2, kappa, type.S)
  #V<- Cii$V
  Psi=varcov.spatial(coords,cov.model=type.S,cov.pars=c(sigmae,phi),nugget=tau2,kappa=kappa)$varcov
  V=(1/sigmae)*Psi

  key <- TRUE
  V.inv <- try(solve(V),silent=TRUE)
  if(class(V.inv) == "try-error") return(Inf)

  detV <- try(determinant(V,logarithm=TRUE)$modulus, silent=TRUE)
  if(class(V.inv) == "try-error") return(Inf)

  AA<- (sum(diag(yy1%*%V.inv))-t(uy1)%*%V.inv%*%media-t(media)%*%V.inv%*%uy1+u1*(t(media)%*%V.inv%*%media))
  vero<- as.numeric(-0.5*(detV+1/sigmae*AA))

  return(-vero)
}

## Funcao para  maximizar usada no SAEM: Estima phi
## tau fixo

FCi.fixo<-function(rhoG,tau2,media,sigmae,kappa,uu,yb,uyb,yyb,type.S,coords){
  phi<-rhoG
  vero=numeric()

  y1 <- yb
  u1 <- uu
  uy1 <- uyb
  yy1 <- yyb

  #Cii<-MatSpatial(MCoord,sigmae, phi, tau2, kappa, type.S)
  #V<- Cii$V
  Psi=varcov.spatial(coords,cov.model=type.S,cov.pars=c(sigmae,phi),nugget=tau2,kappa=kappa)$varcov
  V=(1/sigmae)*Psi

  if(type.S=="exponential"| type.S=="spherical"){
    AA<- (sum(diag(yy1%*%solve(V)))-t(uy1)%*%solve(V)%*%media-t(media)%*%solve(V)%*%uy1+u1*(t(media)%*%solve(V)%*%media))
    vero<- -0.5*(log(det(V))+1/sigmae*AA)
  }

  else{
    key <- TRUE
    V.inv <- try(solve(V),silent=TRUE)
    if(class(V.inv) == "try-error") return(Inf)

    detV <- try(determinant(V,logarithm=TRUE)$modulus, silent=TRUE)
    if(class(V.inv) == "try-error") return(Inf)


    #AA<- (sum(diag(yy1%*%solve(V)))-t(uy1)%*%solve(V)%*%media-t(media)%*%solve(V)%*%uy1+u1*(t(media)%*%solve(V)%*%media))
    AA<- (sum(diag(yy1%*%V.inv))-t(uy1)%*%V.inv%*%media-t(media)%*%V.inv%*%uy1+u1*(t(media)%*%V.inv%*%media))
    vero<- as.numeric(-0.5*(detV+1/sigmae*AA))
  }
  return(-vero)
}


## Log-likelihood ##

LogVerosCens<-function(cc,y,LI,LS,media,Psi,logarithm=FALSE){

  #GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  m<-length(y)

  gammai=media

  if(sum(cc)==0){
    ver<-dmvnorm(as.vector(y),as.vector(gammai),Psi)
  }
  if(sum(cc)>0){
    if(sum(cc)==m){
      auxupper1<-LS-gammai
      auxupper2<-LI-gammai
      ver=prob_opt(lower = c(auxupper2), upper=c(auxupper1), sigma=Psi)
    }
    else{
      muc<- gammai[cc==1,]+Psi[cc==1,cc==0]%*%solve(Psi[cc==0,cc==0])%*%(y[cc==0]-gammai[cc==0,])
      Sc<- Psi[cc==1,cc==1]-Psi[cc==1,cc==0]%*%solve(Psi[cc==0,cc==0])%*%Psi[cc==0,cc==1]
      auxupper1 <- LS[cc==1]-muc
      auxupper2 <- LI[cc==1]-muc
      if(logarithm){ver<-dmvnorm(y[cc==0],gammai[cc==0,],Psi[cc==0,cc==0],log=logarithm)+
        prob_opt(lower=c(auxupper2),upper=c(auxupper1),sigma=Sc,uselog2 = TRUE)/log2(exp(1))} #changing log base from 2 to e
        else{ver<-dmvnorm(y[cc==0],gammai[cc==0,],Psi[cc==0,cc==0])*(prob_opt(lower=c(auxupper2),upper=c(auxupper1),sigma=Sc))}

    }
  }


  obj.out <- list(ver = ver)
  return(obj.out)

}

################################################################################
## Amostrador Gibbs
################################################################################

amostradordegibbs <- function(M,M0,nj,t1,t2,cc1,y1,LI,LS,media,Gama){

  t2 <- 1

  draws <- matrix(NA,nrow=M,ncol=(nj+1))
  draws[1,1:nj] <- t1
  draws[1,nj+1] <- 1

  gammai <- media

  if(sum(cc1)==0){
    for(i in 2:M){
      t1 <- y1
      t2 <- 1
      draws[i,1:nj] <- t1
      draws[i,nj+1] <- t2
    }
  }
  if(sum(cc1)>0 & sum(cc1)==nj){
    for(i in 2:M){
      g <- gammai
      g <- as.vector(g)
      t1 <- as.vector(rtmvnorm(1, mean = g, sigma = ((1/t2)*Gama), lower=LI[cc1==1], upper=LS[cc1==1],
                               algorithm="gibbs", thinning=2))
      t2 <- 1
      draws[i,1:nj] <- t1
      draws[i,nj+1] <- t2
    }
  }
  if(sum(cc1)>0 & sum(cc1)<nj){
    if(sum(cc1)==1){
      for(i in 2:M){
        g <- gammai
        t1[cc1==0] <- y1[cc1==0]
        muc <- g[cc1==1]+Gama[cc1==1,cc1==0]%*%solve(Gama[cc1==0,cc1==0])%*%(y1[cc1==0]-g[cc1==0])
        muc <- as.vector(muc)
        Sc <- Gama[cc1==1,cc1==1]-Gama[cc1==1,cc1==0]%*%solve(Gama[cc1==0,cc1==0])%*%Gama[cc1==0,cc1==1]
        Sc <- as.numeric(Sc)
        y_r <-  rtnorm(1, mean = muc, sd=(sqrt(Sc/t2)),  lower=LI[cc1==1], upper=LS[cc1==1])
        t1[cc1==1] <- y_r
        t2 <- 1
        draws[i,1:nj] <- t1
        draws[i,nj+1] <- t2
      }
    }
    else{
      for(i in 2:M){
        g <- gammai
        t1[cc1==0] <- y1[cc1==0]
        muc <- g[cc1==1]+Gama[cc1==1,cc1==0]%*%solve(Gama[cc1==0,cc1==0])%*%(y1[cc1==0]-g[cc1==0])
        muc <- as.vector(muc)
        Sc <- Gama[cc1==1,cc1==1]-Gama[cc1==1,cc1==0]%*%solve(Gama[cc1==0,cc1==0])%*%Gama[cc1==0,cc1==1]
        y_r <- rtmvnorm(1, mean = muc, sigma = ((1/t2)*Sc), lower=LI[cc1==1], upper=LS[cc1==1],
                        algorithm="gibbs", thinning=2)
        t1[cc1==1] <- y_r
        t2 <- 1
        draws[i,1:nj] <- t1
        draws[i,nj+1] <- t2
      }
    }
  }


  # Amostra com burnin (M0)
  amostragibbs <- draws[(M0+1):M,]

  obj.out <- list(amostragibbs = amostragibbs)
  return(obj.out)

}


#############

dist<-function(coords){

  n=nrow(coords)
  dist=matrix(0,n,n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      dist[i,j] = sqrt(sum((coords[i,]-coords[j,])^2));
    }
  }
  dist=dist+t(dist);
  return(list(dist=dist))
}

###############


derivater=function(coords,phi,kappa,type){
  H=dist(coords)
  H=abs(H$dist)
  if(type=="exponential"){
    H1 <- (H/phi^2)*exp(-(H/phi))
    diag(H1) <- 0

    H2 <- ((H*(H-(2*phi)))/phi^4)*exp(-(H/phi))
    diag(H2) <- 0
  }

  if(type=="gaussian"){
    H1 <- (2/phi^3)*H^2*exp(-((H/phi)^2))
    diag(H1) <- 0
    H2 <- (((4*H^4)-(6*H^2*phi^2))/phi^6)*exp(-((H/phi)^2))
    diag(H1) <- 0
  }

  if(type=="matern"){

    Ak=besselK((H/phi),nu=(kappa-1))+besselK((H/phi),nu=(kappa+1))
    Bk=besselK((H/phi),nu=(kappa-2))-(2*besselK((H/phi),nu=(kappa)))+besselK((H/phi),nu=(kappa+2))
    cons=-1/((2^kappa)*phi^2*gamma(kappa))

    H1=cons*((H/phi)^kappa)*((H*Ak)+(2*kappa*phi*besselK((H/phi),nu=(kappa))))
    diag(H1) <- 0
    t1=-4*phi*(kappa+1)*((kappa*phi*besselK((H/phi),nu=(kappa)))+(H*Ak))
    t2=-(H^2/(4*phi^6))*((H/phi)^(kappa-2))*(-H^2*Bk)
    H2=(2^(1-kappa)/gamma(kappa))*(t1+t2)
    diag(H2) <- 0
  }

  if(type=="spherical"){
    H1=matrix(0,dim(H)[1],dim(H)[1])
    H2=matrix(0,dim(H)[1],dim(H)[1])
    H1[H<phi]=(1.5*H[H<phi]/phi^2)-(1.5*(H[H<phi]^3)/phi^4)
    H2[H<phi]=(-3*H[H<phi]/phi^3)+(6*(H[H<phi]^3)/phi^5)
    diag(H1) <- 0
    diag(H2) <- 0
  }


  if(type=="powered.exponential"){
    H1=(kappa*H^kappa/(phi^(kappa+1)))*exp(-((H/phi)^kappa))
    H2=H1*((kappa*H^kappa/(phi^(kappa+1)))-((kappa+1)/phi))
    diag(H1) <- 0
    diag(H2) <- 0
  }

  if(type=="cauchy"){
    H1=((-2*kappa/(phi^3))*(H^2))*((1+ (H/phi)^2)^(kappa-1))
    te1=((3/(phi^4))*(H^2))*((1+ (H/phi)^2)^(kappa-1))
    te2=(((2*(kappa-1))/(phi^6))*(H^4))*((1+ (H/phi)^2)^(kappa-2))
    H2=(-2*kappa)*(te1-te2)
    diag(H1) <- 0
    diag(H2) <- 0
  }


  devR1 <- H1
  devR2 <- H2
  return(list(d1=devR1,d2=devR2,H=H))
}

################################################

derivQ=function(est,fix.nugget=T){
  qlog=function(theta,coords,X,cov.model,kappa,uy,uyy){
    k=dim(X)[2]
    beta=theta[1:k]
    if(fix.nugget==T){
      tau2=est$tau2
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

  beta=est$beta1
  sigma2=est$sigmae
  phi=est$phi
  tau2=est$tau2
  if(fix.nugget==T){
    theta=c(beta,sigma2,phi)
  }
  else{
    theta=c(beta,sigma2,phi,tau2)
  }

  type=est$type
  uy=est$uy
  uyy=est$yy
  beta=est$beta1
  X=as.matrix(est$X)
  coords=est$coords
  Psi<-varcov.spatial(coords,cov.model=type,cov.pars=c(sigma2,phi),nugget=tau2,kappa=kappa)$varcov
  r=qlog(theta,coords=coords,X=X,cov.model=type,kappa=kappa,uy=uy,uyy=uyy)
  s=grad(qlog,theta,coords=coords,X=X,cov.model="matern",kappa=kappa,uy=uy,uyy=uyy)
  Q=hessian(qlog,theta,coords=coords,X=X,cov.model="matern",kappa=kappa,uy=uy,uyy=uyy)
  Qinv=solve(-Q)

  return(list(Qlogvalue=r,gradQ=s,HQ=Q,QI=Qinv,Sigma=Psi))

}


###########################################################


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
  tau2=est$tau2
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



