###########
##########2.Seminaive algortihm for spatial censored prediction############
###########

seminaive=function(data,y.col,coords.col,covar,covar.col,copred,cov.model,thetaini,fix.nugget=T,nugget,kappa,cons,MaxIter,cc,cutof,trend){
  if(covar==T){
    geodataprin=as.geodata(data,coords.col =coords.col, data.col = y.col, covar.col = covar.col)
    covariates=geodataprin$covariate
    covparini=thetaini[1:2]
  }

  if(covar==F){
    geodataprin=as.geodata(data,coords.col =coords.col, data.col = y.col)
    covparini=thetaini
  }


  y=geodataprin$data
  #X=as.matrix(rep(1,length(y)))
  zk=y
  zk[cc==1]=0

  coords=geodataprin$coords
  sigmau=var(y)
  skewu=skewness(y)
  crit1=crit2=crit3=1
  count=0

  while(crit1>cons[1] | crit2>(cons[2]*sigmau) | crit3<(cons[3]*skewu)){
    count=count+1
    if(covar==T){
      dados=cbind(coords,zk,covariates)
      geodata=as.geodata(dados,coords.col=coords.col,y.col=y.col,covar.col=covar.col)
    }

    if(covar==F){
      dados=cbind(coords,zk)
      geodata=as.geodata(dados,coords.col=coords.col,y.col=y.col)
    }
    #variog=variog(geodata)
    #s=variofit(variog,ini.cov.pars=covparini, cov.model=cov.model,fix.nugget =fix.nugget, nugget =nugget,fix.kappa = TRUE, kappa = kappa,weights="cressie")

    s=likfit(geodata, coords = geodata$coords, data = geodata$data,
             trend = trend, ini.cov.pars=covparini, fix.nugget = fix.nugget, nugget = nugget,cov.model=cov.model,kappa=kappa)

    phi=s$cov.pars[2]
    sigma2=s$cov.pars[1]
    tau2=s$nugget
    beta=s$beta
    AIC=s$AIC
    BIC=s$BIC
    loglik=s$loglik
    #V=varcov.spatial(geodata$coords,cov.model="matern",cov.pars=c(sigma2,phi),nugget=tau2,kappa=0.5)$varcov
    theta=c(beta,sigma2,phi,tau2)

    z=0
    v=0

    for(i in 1:sum(cc==1)){
      coord=geodata$coords[-i,]
      ysin=geodata$data[-i]
      if(covar==T){
        k=dim(covariates)[2]
        covsin=geodata$covariate[-i,1:2]
        u=as.geodata(cbind(coord,ysin,covsin),coords.col=coords.col,y.col=y.col,covar.col=covar.col)
      }

      if(covar==F){
        k=1
        u=as.geodata(cbind(coord,ysin),coords.col=coords.col,y.col=y.col)

      }

      kc=krige.control(type.krige = "ok",cov.model=cov.model, obj.model=s, kappa=kappa)

      a=krige.conv(u,coords=coord, data=ysin,locations=geodata$coords[i,], krige=kc,output=output.control(messages=F))$predict
      z[i]=krige.conv(u,coords=coord, data=ysin,locations=geodata$coords[i,], krige=kc,output=output.control(messages=F))$predict
      v[i]=min(cutof[i],z[i])
    }
    zk1=geodata$data
    zk1[cc==1]=v
    if (count>0){
      crit1=abs((var(zk1)-var(zk))/var(zk))
      crit2=var(zk1)
      crit3=skewness(zk1)
    }
    if (count==MaxIter){
      crit1 <-0.000000000000001
      crit2 <- 0.000000000000001
      crit3 <- 100000

    }
    zk=zk1
  }

  kc=krige.control(type.krige = "ok",obj.model=s, kappa=kappa,
                   nugget=theta[(k+3)])
  ycons=zk
  if(covar==T){
    ycons=data.frame(coords,ycons,covariates)
    ycons=as.geodata(ycons)
  }
  if(covar==F){
    ycons=data.frame(coords,ycons)
    ycons=as.geodata(ycons)
  }

  aux=krige.conv(ycons,coords=coords, data=ycons$data,locations=copred, krige=kc,output=output.control(messages=F))
  predictions=aux$predict
  sdpred=sqrt(aux$krige.var)
  return(list(zk=zk,AIC=AIC,BIC=BIC,beta=beta,theta=theta,predictions=predictions,sdpred=sdpred,loglik=loglik))

}
