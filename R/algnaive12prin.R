###########
####### 3.Naive 1 and Naive 2 method for spatial prediction #####################
###########



naive12=function(data,cc,copred,thetaini,y.col=3,coords.col=1:2,covar,covar.col,fix.nugget=T,nugget,kappa=0.5,cutof,cov.model,trend){

  if(covar==T){
    geodataprin=as.geodata(data,coords.col =coords.col, data.col = y.col, covar.col = covar.col)
  }

  if(covar==F){
    geodataprin=as.geodata(data,coords.col =coords.col, data.col = y.col)
  }

  coords=geodataprin$coords
  y=geodataprin$data
  covparini=thetaini
  naive1=y
  naive2=y
  naive1[cc==1]=cutof
  naive2[cc==1]=cutof/2

  if(covar==T){
    covariates=geodataprin$covariate
    dataa=data.frame(coords,naive1,covariates)
    datab=data.frame(coords,naive2,covariates)
    geodata1=as.geodata(dataa,coords.col =coords.col, data.col = y.col, covar.col = covar.col)
    geodata2=as.geodata(datab,coords.col =coords.col, data.col = y.col, covar.col = covar.col)

  }

  if(covar==F){
    dataa=data.frame(coords,naive1)
    datab=data.frame(coords,naive2)
    geodata1=as.geodata(dataa,coords.col =coords.col, data.col = y.col)
    geodata2=as.geodata(datab,coords.col =coords.col, data.col = y.col)
  }



  est1=likfit(geodata1, coords = geodata1$coords, data = geodata1$data,kappa=kappa,
              trend = trend, ini.cov.pars=covparini, fix.nugget = fix.nugget, nugget = nugget,cov.model=cov.model)

  est2=likfit(geodata2, coords = geodata2$coords, data = geodata2$data,kappa=kappa,
              trend =  trend, ini.cov.pars=covparini, fix.nugget = fix.nugget, nugget =nugget,cov.model=cov.model)
  kc1=krige.control(type.krige = "ok",obj.model=est1)
  kc2=krige.control(type.krige = "ok",obj.model=est2)

  phi1=est1$cov.pars[2]
  sigma21=est1$cov.pars[1]
  tau21=est1$nugget
  beta1=est1$beta
  #V=varcov.spatial(geodata$coords,cov.model="matern",cov.pars=c(sigma2,phi),nugget=tau2,kappa=0.5)$varcov
  theta1=c(beta1,sigma21,phi1,tau21)

  phi2=est2$cov.pars[2]
  sigma22=est2$cov.pars[1]
  tau22=est2$nugget
  beta2=est2$beta
  #V=varcov.spatial(geodata$coords,cov.model="matern",cov.pars=c(sigma2,phi),nugget=tau2,kappa=0.5)$varcov
  theta2=c(beta2,sigma22,phi2,tau22)
  AIC1=est1$AIC
  BIC1=est1$BIC
  AIC2=est2$AIC
  BIC2=est2$BIC
  loglik1=est1$loglik
  loglik2=est2$loglik
  aux1=krige.conv(geodata1,coords=coords, data=geodata1$data,locations=copred, krige=kc1,output=output.control(messages=F))
  aux2=krige.conv(geodata2,coords=coords, data=geodata2$data,locations=copred, krige=kc2,output=output.control(messages=F))
  predictions1=aux1$predict
  predictions2=aux2$predict
  sdpred1=sqrt(aux1$krige.var)
  sdpred2=sqrt(aux2$krige.var)

  return(list(beta1=beta1,beta2=beta2,theta1=theta1,theta2=theta2,predictions1=predictions1,predictions2=predictions2,sdpred1=sdpred1,sdpred2=sdpred2,AIC1=AIC1,BIC1=BIC1,AIC2=AIC2,BIC2=BIC2,loglik1=loglik1,loglik2=loglik2))


}
