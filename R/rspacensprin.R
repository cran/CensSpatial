##############
#############1.simcenspadata:  Censored Spatial data simulation#######################################
##############

simcenspadata=function(cov.pars,beta,x=as.matrix(rep(1,n)),coords,kappa,cens,n,n1,type,cens.type){
  phi=cov.pars[1]
  sigma2=cov.pars[2]
  tau2=cov.pars[3]
  cc<-rep(0,n-n1)
  x=as.matrix(x)
  if((dim(x)[2])==1){
    media=x*beta
  }else{
    media=x%*%beta
  }
  covar<-varcov.spatial(coords,cov.model=type,cov.pars=c(sigma2,phi),nugget=tau2,kappa=0.5)$varcov

  media1=media
  y<-c(rmvnorm(1,media1,covar))
  nnopred=n-n1

  coords1=coords[(nnopred+1):n,]

  coords2=coords[1:nnopred,]

  y1 <-y[1:nnopred]
  aa =  sort(y1)

  if(cens.type=="left"){
    bb =  aa[1:(cens*(nnopred))]
    cutoff <- bb[(cens*(nnopred))]
    cc =  matrix(1,(n-n1),1)*(y1< cutoff)
  }
  if(cens.type=="right"){
    bb =  aa[1:((1-cens)*(nnopred))]
    cutoff <- bb[((1-cens)*(nnopred))]
    cc =  matrix(1,(n-n1),1)*(y1>cutoff)
  }
  y1[cc==1]=  cutoff
  valre=y[(nnopred+1):n]
  datare=data.frame(coords2,y1)

  s=list(y=y,datare=datare,valre=valre,cc=cc,cutoff=cutoff,coords1=coords1)

  return(s)
}

