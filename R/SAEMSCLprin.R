#######################
##################6.SAEM Algorithm estimation##################
#######################

SAEM_Spatial<-function(cc,y,cens.type,trend,LI,LS,x,coords,kappa,M=20,perc=0.25,MaxIter=300,pc=0.2,type.S,
                       fix.nugget=TRUE,nugget,inits.sigmae,inits.phi,search=F,lower,upper){
  pb <- txtProgressBar(min = 0, max = MaxIter, style = 3)
  #valores iniciais
  m<-length(y)



  if(trend=="cte"){
    x=as.matrix(rep(1,m))
  }

  if(trend=="1st"){
    x=as.matrix(cbind(1,coords))
  }

  interaction=coords[,1]*coords[,2]
  if(trend=="2nd"){
    x=as.matrix(cbind(1,coords,(coords)^2,interaction))
  }

  if(trend=="other"){
    x=as.matrix(x)
  }



  beta1<-solve(t(x)%*%x)%*%t(x)%*%y

  p<-length(beta1)


  if(cens.type=="left"){
    LI=rep(-Inf,length(cc))
    LS=rep(Inf,length(cc))
    LS[cc==1]=y[cc==1]
    LI=as.vector(LI)
    LS=as.vector(LS)
  }

  if(cens.type=="right"){
    LI=rep(-Inf,length(cc))
    LI[cc==1]=y[cc==1]
    LS=rep(Inf,length(cc))
    LI=as.vector(LI)
    LS=as.vector(LS)
  }

  if(cens.type=="both"){
    LI=LI
    LS=LS
    LI=as.vector(LI)
    LS=as.vector(LS)
  }

  # Parametro Tau Fixo
  if(fix.nugget==TRUE){

    tau2<-nugget
    sigmae<-inits.sigmae
    phi<-inits.phi

    teta <- c(beta1,sigmae,phi)

    MG <- round(M/(1-perc),0)
    M0 <- MG - M

    Theta <- matrix(NA,MaxIter,length(teta))

    tyi <- NULL
    tui <- runif(1,0,1) # valor inicial de ui
    tyy <- NULL

    criterio <- 1
    count <- 0

    media= x%*%beta1
    #AuxM<-MatSpatial(MCoord, sigmae, phi, tau2, kappa, type.S)
    #V<-AuxM$V  ### (variancia-Covariancia)/sigmae
    #Psi<-AuxM$Cov   ## variancia-Covariancia
    Psi=varcov.spatial(coords,cov.model=type.S,cov.pars=c(sigmae,phi),nugget=tau2,kappa=kappa)$varcov
    V=(1/sigmae)*Psi

    if(pc==1)
    {
      seqq=rep(1,pc*MaxIter)
    } else
    {
      seqq = c(rep(1,pc*MaxIter),(1/((((pc*MaxIter)+1):MaxIter)-(pc*MaxIter))))
      seqq = c(rep(1,MaxIter-length(seqq)),seqq)
    }

    SAEM_ss <- array(data=0,dim=c(MaxIter+1))
    SAEM_xx <- array(data=0,dim=c(MaxIter+1,p,p))
    SAEM_xy <- array(data=0,dim=c(MaxIter+1,p))

    SAEM_y <- array(data=0,dim=c(MaxIter+1,m))
    SAEM_u <- array(data=0,dim=c(MaxIter+1))
    SAEM_uy <- array(data=0,dim=c(MaxIter+1,m))
    SAEM_yy <- array(data=0,dim=c(MaxIter+1,m,m))

    initime=Sys.time()
    while(criterio > 0.00001){

      count <- count + 1
      # print(count)
      setTxtProgressBar(pb, count)
      tu1 <- tui

      ## Passo de Simulacao: Gera das distribuicoes condicionis
      t1 <- y
      t2 <- tu1
      gibbs <- amostradordegibbs(MG,M0,m,t1,t2,cc,y,LI,LS,media,Psi)
      amostragibbs <- gibbs$amostragibbs

      uyi <- matrix(amostragibbs[,1:m],nrow=M,ncol=m)
      uui <- amostragibbs[,m+1]

      ## Passo de Aproximacao
      somass <- 0
      somaxx <- matrix(0,p,p)
      somaxy <- matrix(0,p,1)
      somay <- matrix(0,m,1)
      somau <- 0
      somauy <- matrix(0,m,1)
      somayy <- matrix(0,m,m)
      for(k in 1:M){
        yi <- matrix(uyi[k,],nrow=m,ncol=1)
        ui <- uui[k]
        somass <- somass + ui*(t(yi-media)%*%solve(V)%*%(yi-media))
        somaxx <- somaxx + ui*(t(x)%*%solve(V)%*%x)
        somaxy <- somaxy + ui*t(x)%*%solve(V)%*%(yi)
        somay <- somay + yi
        somau <- somau + ui
        somauy <- somauy + ui*yi
        somayy <- somayy + ui*yi%*%t(yi)
      }

      E_ss = (1/M)*somass
      E_xx = (1/M)*somaxx
      E_xy = (1/M)*somaxy
      E_y = (1/M)*somay
      E_u = (1/M)*somau
      E_uy = (1/M)*somauy
      E_yy = (1/M)*somayy

      ## Aproximacao Estocastica
      SAEM_ss[count+1] = SAEM_ss[count] + seqq[count]*(E_ss - SAEM_ss[count])
      SAEM_xx[count+1,,] = SAEM_xx[count,,] + seqq[count]*(E_xx - SAEM_xx[count,,])
      SAEM_xy[count+1,] = SAEM_xy[count,] + seqq[count]*(E_xy - SAEM_xy[count,])

      SAEM_y[count+1,] = SAEM_y[count,] + seqq[count]*(E_y - SAEM_y[count,])
      SAEM_u[count+1] = SAEM_u[count] + seqq[count]*(E_u - SAEM_u[count])
      SAEM_uy[count+1,] = SAEM_uy[count,] + seqq[count]*(E_uy - SAEM_uy[count,])
      SAEM_yy[count+1,,] = SAEM_yy[count,,] + seqq[count]*(E_yy - SAEM_yy[count,,])

      tyi<- SAEM_y[count+1,]
      tui<-SAEM_u[count+1]
      tyy<-SAEM_yy[count+1,,]

      ## Passo M
      beta1<- solve(SAEM_xx[count+1,,])%*%SAEM_xy[count+1,]
      media<-x%*%beta1

      sigmae<- (1/m)*(SAEM_ss[count+1])
      sigmae<-as.numeric(sigmae)

      if(search==F){
        lower=0
        upper=10000
      }

      if(search==T){
        lower=lower
        upper=upper
      }

      phi <- optimize(f=FCi.fixo, lower = lower, upper = upper, tau2 = tau2, media = media, sigmae = sigmae,
                      kappa = kappa, uu = SAEM_u[count+1], yb = SAEM_y[count+1,], uyb = SAEM_uy[count+1,],
                      yyb = SAEM_yy[count+1,,], type.S = type.S,coords = coords)$minimum
      teta1 <- c(beta1,sigmae,phi)
      #  print(teta1)

      # AuxM<-MatSpatial(MCoord, sigmae, phi ,tau2 , kappa, type.S)
      #V<-AuxM$V
      #Psi<-AuxM$Cov
      Psi=varcov.spatial(coords,cov.model=type.S,cov.pars=c(sigmae,phi),nugget=tau2,kappa=kappa)$varcov
      V=(1/sigmae)*Psi

      if (count>1){
        criterio <- sqrt((teta1/teta-1)%*%(teta1/teta-1))
      }

      if (count==MaxIter){
        criterio <- 0.00000000001
      }

      Theta[count,] <- teta1
      teta<-teta1

    }
    endtime=Sys.time()
    timediffe=endtime-initime
    Theta=Theta[1:count,]
    logver <- log(LogVerosCens(cc,y,LI,LS,media,Psi)$ver)
    npar<-length(c(teta1))
    loglik<-logver
    AICc<- -2*loglik +2*npar
    AICcorr<- AICc + ((2*npar*(npar+1))/(m-npar-1))
    BICc <- -2*loglik +log(m)*npar

  }

  # Parametro Tau Estimado
  if(fix.nugget==FALSE){

    tau2<-nugget
    sigmae<-inits.sigmae
    phi<-inits.phi

    teta <- c(beta1,sigmae,phi,tau2)

    MG <- round(M/(1-perc),0)
    M0 <- MG - M

    Theta <- matrix(NA,MaxIter,length(teta))

    tyi <- NULL
    tui <- runif(1,0,1) # valor inicial de ui
    tyy <- NULL

    criterio <- 1
    count <- 0

    media= x%*%beta1
    #AuxM<-MatSpatial(MCoord, sigmae, phi, tau2, kappa, type.S)
    #V<-AuxM$V  ### (variancia-Covariancia)/sigmae
    #Psi<-AuxM$Cov   ## variancia-Covariancia
    Psi=varcov.spatial(coords,cov.model=type.S,cov.pars=c(sigmae,phi),nugget=tau2,kappa=kappa)$varcov
    V=(1/sigmae)*Psi
    if(pc==1)
    {
      seqq=rep(1,pc*MaxIter)
    } else
    {
      seqq = c(rep(1,pc*MaxIter),(1/((((pc*MaxIter)+1):MaxIter)-(pc*MaxIter))))
      seqq = c(rep(1,MaxIter-length(seqq)),seqq)
    }

    SAEM_ss <- array(data=0,dim=c(MaxIter+1))
    SAEM_xx <- array(data=0,dim=c(MaxIter+1,p,p))
    SAEM_xy <- array(data=0,dim=c(MaxIter+1,p))

    SAEM_y <- array(data=0,dim=c(MaxIter+1,m))
    SAEM_u <- array(data=0,dim=c(MaxIter+1))
    SAEM_uy <- array(data=0,dim=c(MaxIter+1,m))
    SAEM_yy <- array(data=0,dim=c(MaxIter+1,m,m))

    initime=Sys.time()

    while(criterio > 0.00001){
      setTxtProgressBar(pb, count)
      count <- count + 1
      #print(count)

      tu1 <- tui

      ## Passo de Simulacao: Gera das distribuicoes condicionis
      t1 <- y
      t2 <- tu1
      gibbs <- amostradordegibbs(MG,M0,m,t1,t2,cc,y,LI,LS,media,Psi)
      amostragibbs <- gibbs$amostragibbs

      uyi <- matrix(amostragibbs[,1:m],nrow=M,ncol=m)
      uui <- amostragibbs[,m+1]

      ## Passo de Aproximacao
      somass <- 0
      somaxx <- matrix(0,p,p)
      somaxy <- matrix(0,p,1)
      somay <- matrix(0,m,1)
      somau <- 0
      somauy <- matrix(0,m,1)
      somayy <- matrix(0,m,m)
      for(k in 1:M){
        yi <- matrix(uyi[k,],nrow=m,ncol=1)
        ui <- uui[k]
        somass <- somass + ui*(t(yi-media)%*%solve(V)%*%(yi-media))
        somaxx <- somaxx + ui*(t(x)%*%solve(V)%*%x)
        somaxy <- somaxy + ui*t(x)%*%solve(V)%*%(yi)
        somay <- somay + yi
        somau <- somau + ui
        somauy <- somauy + ui*yi
        somayy <- somayy + ui*yi%*%t(yi)
      }

      E_ss = (1/M)*somass
      E_xx = (1/M)*somaxx
      E_xy = (1/M)*somaxy
      E_y = (1/M)*somay
      E_u = (1/M)*somau
      E_uy = (1/M)*somauy
      E_yy = (1/M)*somayy

      ## Aproximacao Estocastica
      SAEM_ss[count+1] = SAEM_ss[count] + seqq[count]*(E_ss - SAEM_ss[count])
      SAEM_xx[count+1,,] = SAEM_xx[count,,] + seqq[count]*(E_xx - SAEM_xx[count,,])
      SAEM_xy[count+1,] = SAEM_xy[count,] + seqq[count]*(E_xy - SAEM_xy[count,])

      SAEM_y[count+1,] = SAEM_y[count,] + seqq[count]*(E_y - SAEM_y[count,])
      SAEM_u[count+1] = SAEM_u[count] + seqq[count]*(E_u - SAEM_u[count])
      SAEM_uy[count+1,] = SAEM_uy[count,] + seqq[count]*(E_uy - SAEM_uy[count,])
      SAEM_yy[count+1,,] = SAEM_yy[count,,] + seqq[count]*(E_yy - SAEM_yy[count,,])

      tyi<- SAEM_y[count+1,]
      tui<-SAEM_u[count+1]
      tyy<-SAEM_yy[count+1,,]

      ## Passo M
      beta1<- solve(SAEM_xx[count+1,,])%*%SAEM_xy[count+1,]
      media<-x%*%beta1

      sigmae<- (1/m)*(SAEM_ss[count+1])
      sigmae<-as.numeric(sigmae)

      if(search==F){
        lower=lower=c(0,0)
        upper=c(10000,10000)
      }

      if(search==T){
        lower=lower
        upper=upper
      }


      #rhos <- optim(c(phi,tau2), method = "L-BFGS-B", FCi, lower =, upper =upper, media = media, sigmae = sigmae,
      #              MCoord = MCoord, kappa = kappa, uu = SAEM_u[count+1], yb = SAEM_y[count+1,], uyb = SAEM_uy[count+1,],
      #              yyb = SAEM_yy[count+1,,], type.S = type.S, hessian = TRUE, control=c(trace=5))$par

      rhos <- optimx(c(phi,tau2), method = "nlminb", fn=FCi, lower=lower, upper=upper, media = media, sigmae = sigmae,
                     kappa = kappa, uu = SAEM_u[count+1], yb = SAEM_y[count+1,], uyb = SAEM_uy[count+1,],
                     yyb = SAEM_yy[count+1,,], type.S = type.S, coords = coords,hessian = T)
      phi<-rhos$p1
      tau2<-rhos$p2
      teta1 <- c(beta1,sigmae,phi,tau2)
      # print(teta1)

      #AuxM<-MatSpatial(MCoord, sigmae, phi ,tau2 , kappa, type.S)
      #V<-AuxM$V
      #  Psi<-AuxM$Cov
      Psi=varcov.spatial(coords,cov.model=type.S,cov.pars=c(sigmae,phi),nugget=tau2,kappa=kappa)$varcov
      V=(1/sigmae)*Psi
      if (count>1){
        criterio <- sqrt((teta1/teta-1)%*%(teta1/teta-1))
      }

      if (count==MaxIter){
        criterio <- 0.00000000001
      }

      Theta[count,] <- teta1
      teta<-teta1
    }

    endtime=Sys.time()
    Theta=Theta[1:count,]

    logver <- LogVerosCens(cc,y,LI,LS,media,Psi,logarithm=TRUE)$ver
    npar<-length(c(teta1))
    loglik<-logver
    AICc<- -2*loglik +2*npar
    AICcorr<- AICc + ((2*npar*(npar+1))/(m-npar-1))
    BICc <- -2*loglik +log(m)*npar
    setTxtProgressBar(pb, MaxIter)
    timediffe=endtime-initime
  }

  fitted=x%*%beta1
  obj.out <- list(beta1 = beta1, sigmae = sigmae, phi = phi, tau2 = tau2, Theta=Theta, loglik=loglik,
                  AIC=AICc, BIC=BICc, AICcorr=AICcorr,X=x, Psi=Psi,trend=trend,
                  theta = teta1, yy = tyy,uy=tyi,cc=cc,type=type.S,kappa=kappa,coords=coords,timex=timediffe,iter=count,fitted=fitted)

  class(obj.out) <- "SAEMSpatialCens"

  return(obj.out)

}
