
###########10.Local influence measures#######


locinme=function(est,fix.nugget,diag.plot=TRUE,type.plot="all",c=3){

  beta=est$beta
  coords = est$coords
  type = est$type
  sigma2=est$sigma2
  phi=est$phi
  theta=c(beta,sigma2,phi)
  uy=est$uy
  uyy=est$uyy
  beta=est$beta
  n=length(uy)
  k=length(beta)
  kappa=est$kappa
  pb <- txtProgressBar(min = 0, max =n, style = 3)
  X=as.matrix(est$X)
  tau2=est$nugget
  Sigma=est$Psi
  V1=solve(Sigma)
  V=(1/sigma2)*Sigma
  rphi=V-((tau2/sigma2)*diag(1,dim(V)[1],dim(V)[1]))


  if(fix.nugget==T){
    d1=derivQ(est)
    derive=derivater(coords,phi=phi,kappa=kappa,type=type)
    deri1phi=sigma2*derive$d1
    deri1sigma=Sigma
    iden=diag(1,n)
    P=X%*%solve(t(X)%*%V1%*%X)%*%t(X)%*%V1
    deri1invphi=-V1%*%deri1phi%*%V1

    deri1invsigma=-V1%*%rphi%*%V1

    ##response perturbation#####

    deltabeta1=-t(X)%*%Sigma

    deltaphi1=0
    deltasigma1=0
    multiphi=t(P%*%uy)%*%deri1invphi
    multisigma=t(P%*%uy)%*%deri1invsigma

    ####scale matrix perturbation############
    deltabeta2=matrix(0,k,n)
    deltaalpha2=matrix(0,2,n)
    D=iden

    ######explanatory variable perturbation######

    deltabeta3=matrix(0,n,k)
    deltaalpha3=matrix(0,2,n)
    betacal=solve(t(X)%*%V1%*%X)%*%t(X)%*%V1%*%uy


    for(i in 1:n){

      setTxtProgressBar(pb,i)

      ##response perturbation#####

      Z=matrix(0,n,n)
      Z[,i]=uy
      deltaphi1[i]=(0.5*tr((t(Z)+Z)%*%deri1invphi))+ multiphi[i]
      deltasigma1[i]=(0.5*tr((t(Z)+Z)%*%deri1invsigma))+ multisigma[i]


      ####scale matrix perturbation############

      d=matrix(0,n,n)
      d[i,i]=1
      deltabeta2[,i]=t(X)%*%V1%*%d%*%(iden-P)%*%uy

      a1=-0.5*tr((V1%*%d%*%D%*%deri1phi)+(V1%*%D%*%d%*%deri1phi))
      b11=tr(uyy%*%V1%*%deri1phi%*%V1%*%d)
      b12=t(uy)%*%(t(P)-(2*iden))%*%(V1%*%deri1phi%*%V1)%*%(d%*%P%*%uy)
      b1=-0.5*(b11+b12)
      deltaphi2=a1+b1


      c1=-0.5*tr((V1%*%d%*%D%*%deri1sigma)+(V1%*%D%*%d%*%deri1sigma))
      e11=tr(uyy%*%V1%*%deri1sigma%*%V1%*%d)
      e12=t(uy)%*%(t(P)-(2*iden))%*%(V1%*%deri1sigma%*%V1)%*%(d%*%P%*%uy)
      e1=-0.5*(e11+e12)
      deltasigma2=c1+e1

      deltaalpha2[,i]=rbind(deltaphi2,deltasigma2)


      ######explanatory variable perturbation######

      wj=matrix(0,n,k)
      wj[i,]=1
      deltabeta3[i,]=t(uy)%*%(iden-(2*P))%*%V1%*%wj
      deltaphi3=t(uy)%*%(iden-t(P))%*%deri1invphi%*%wj%*%betacal
      deltasigma3=t(uy)%*%(iden-t(P))%*%deri1invsigma%*%wj%*%betacal
      deltaalpha3[,i]=c(deltaphi3,deltasigma3)


    }

    ##response perturbation#####

    deltaomega1=rbind(deltabeta1,deltaphi1,deltasigma1)

    Qw1=t(deltaomega1)%*%d1$QI%*%deltaomega1
    desc1=eigen(2*Qw1)
    desc1$values[desc1$values<1e-20]=0
    erara1=desc1$values/sum(desc1$values)
    m01=t(erara1*t((desc1$vectors)^2))
    m01=apply(m01,1,FUN=sum)

    #Bfqd=0
    #for(i in 1:n){
    #dl=rep(0,200)
    #dl[i]=1
    #Bfqd[i]=2*t(dl)%*%Qw1%*%dl/tr(2*Qw1)
    #}

    #lim=mean(m01)+(3*sd(m01))

    ####scale matrix perturbation############

    deltaomega2=rbind(deltabeta2,deltaalpha2)

    Qw2=t(deltaomega2)%*%d1$QI%*%deltaomega2
    desc2=eigen(2*Qw2)
    desc2$values[desc2$values<1e-20]=0
    erara2=desc2$values/sum(desc2$values)
    m02=t(erara2*t((desc2$vectors)^2))
    m02=apply(m02,1,FUN=sum)

    #Bfqd=0
    #for(i in 1:n){
    #dl=rep(0,200)
    #dl[i]=1
    #Bfqd[i]=2*t(dl)%*%Qw%*%dl/tr(2*Qw)
    #}

    #lim=mean(m0)+(3*sd(m0))



    ######explanatory variable perturbation######

    deltaomega3=rbind(t(deltabeta3),deltaalpha3)

    Qw3=t(deltaomega3)%*%d1$QI%*%deltaomega3
    desc3=eigen(2*Qw3)
    desc3$values[desc3$values<1e-20]=0
    erara3=desc3$values/sum(desc3$values)
    m03=t(erara3*t((desc3$vectors)^2))
    m03=apply(m03,1,FUN=sum)

    #Bfqd=0
    #for(i in 1:n){
    #dl=rep(0,200)
    #dl[i]=1
    #Bfqd[i]=2*t(dl)%*%Qw%*%dl/tr(2*Qw)
    #}

    #lim=mean(m0)+(3*sd(m0))



  }

  if(fix.nugget==F){
    d1=derivQ(est,fix.nugget=F)
    derive=derivater(coords,phi=phi,kappa=kappa,type=type)
    deri1phi=sigma2*derive$d1
    deri1sigma=Sigma
    deri1tau=sigma2*diag(1,n)
    iden=diag(1,n)
    P=X%*%solve(t(X)%*%V1%*%X)%*%t(X)%*%V1
    deri1invphi=-V1%*%deri1phi%*%V1

    deri1invsigma=-V1%*%rphi%*%V1

    deri1invtau=-sigma2*V1%*%V1

    ##response perturbation#####
    deltabeta1=-t(X)%*%Sigma
    deltaphi1=0
    deltasigma1=0
    deltatau1=0
    multiphi=t(P%*%uy)%*%deri1invphi
    multisigma=t(P%*%uy)%*%deri1invsigma
    multitau=t(P%*%uy)%*%deri1invtau

    ####scale matrix perturbation############
    deltabeta2=matrix(0,k,n)
    deltaalpha2=matrix(0,3,n)
    D=iden

    ######explanatory variable perturbation######

    deltabeta3=matrix(0,n,k)
    deltaalpha3=matrix(0,3,n)
    betacal=solve(t(X)%*%V1%*%X)%*%t(X)%*%V1%*%uy

    for(i in 1:n){

      setTxtProgressBar(pb,i)

      ##response perturbation#####

      Z=matrix(0,n,n)
      Z[,i]=uy
      deltaphi1[i]=(0.5*tr((t(Z)+Z)%*%deri1invphi))+ multiphi[i]
      deltasigma1[i]=(0.5*tr((t(Z)+Z)%*%deri1invsigma))+ multisigma[i]
      deltatau1[i]=(0.5*tr((t(Z)+Z)%*%deri1invtau))+ multitau[i]

      ####scale matrix perturbation############

      d=matrix(0,n,n)
      d[i,i]=1
      deltabeta2[,i]=t(X)%*%V1%*%d%*%(iden-P)%*%uy

      a1=-0.5*tr((V1%*%d%*%D%*%deri1phi)+(V1%*%D%*%d%*%deri1phi))
      b11=tr(uyy%*%V1%*%deri1phi%*%V1%*%d)
      b12=t(uy)%*%(t(P)-(2*iden))%*%(V1%*%deri1phi%*%V1)%*%(d%*%P%*%uy)
      b1=-0.5*(b11+b12)
      deltaphi2=a1+b1


      c1=-0.5*tr((V1%*%d%*%D%*%deri1sigma)+(V1%*%D%*%d%*%deri1sigma))
      e11=tr(uyy%*%V1%*%deri1sigma%*%V1%*%d)
      e12=t(uy)%*%(t(P)-(2*iden))%*%(V1%*%deri1sigma%*%V1)%*%(d%*%P%*%uy)
      e1=-0.5*(e11+e12)
      deltasigma2=c1+e1

      f1=-0.5*tr((V1%*%d%*%D%*%deri1tau)+(V1%*%D%*%d%*%deri1tau))
      g11=tr(uyy%*%V1%*%deri1tau%*%V1%*%d)
      g12=t(uy)%*%(t(P)-(2*iden))%*%(V1%*%deri1tau%*%V1)%*%(d%*%P%*%uy)
      g1=-0.5*(g11+g12)
      deltatau2=f1+g1

      deltaalpha2[,i]=rbind(deltaphi2,deltasigma2,deltatau2)

      ######explanatory variable perturbation######

      wj=matrix(0,n,k)
      wj[i,]=1
      deltabeta3[i,]=t(uy)%*%(iden-(2*P))%*%V1%*%wj
      deltaphi3=t(uy)%*%(iden-t(P))%*%deri1invphi%*%wj%*%betacal
      deltasigma3=t(uy)%*%(iden-t(P))%*%deri1invsigma%*%wj%*%betacal
      deltatau3=t(uy)%*%(iden-t(P))%*%deri1invtau%*%wj%*%betacal
      deltaalpha3[,i]=c(deltaphi3,deltasigma3,deltatau3)

    }

    ##response perturbation#####

    deltaomega1=rbind(deltabeta1,deltaphi1,deltasigma1,deltatau1)

    Qw1=t(deltaomega1)%*%d1$QI%*%deltaomega1
    desc1=eigen(2*Qw1)
    desc1$values[desc1$values<1e-20]=0

    erara1=desc1$values/sum(desc1$values)
    m01=t(erara1*t((desc1$vectors)^2))
    m01=apply(m01,1,FUN=sum)

    #Bfqd=0
    #for(i in 1:n){
    #dl=rep(0,200)
    #dl[i]=1
    #Bfqd[i]=2*t(dl)%*%Qw1%*%dl/tr(2*Qw1)
    #}

    #lim=mean(m01)+(3*sd(m01))


    ####scale matrix perturbation############

    deltaomega2=rbind(deltabeta2,deltaalpha2)

    Qw2=t(deltaomega2)%*%d1$QI%*%deltaomega2
    desc2=eigen(2*Qw2)
    desc2$values[desc2$values<1e-20]=0
    erara2=desc2$values/sum(desc2$values)
    m02=t(erara2*t((desc2$vectors)^2))
    m02=apply(m02,1,FUN=sum)

    #Bfqd=0
    #for(i in 1:n){
    #dl=rep(0,200)
    #dl[i]=1
    #Bfqd[i]=2*t(dl)%*%Qw%*%dl/tr(2*Qw)
    #}

    #lim=mean(m0)+(3*sd(m0))

    ######explanatory variable perturbation######

    deltaomega3=rbind(t(deltabeta3),deltaalpha3)

    Qw3=t(deltaomega3)%*%d1$QI%*%deltaomega3
    desc3=eigen(2*Qw3)
    desc3$values[desc3$values<1e-20]=0
    erara3=desc3$values/sum(desc3$values)
    m03=t(erara3*t((desc3$vectors)^2))
    m03=apply(m03,1,FUN=sum)

    #Bfqd=0
    #for(i in 1:n){
    #dl=rep(0,200)
    #dl[i]=1
    #Bfqd[i]=2*t(dl)%*%Qw%*%dl/tr(2*Qw)
    #}

    #lim=mean(m0)+(3*sd(m0))
  }

  lim1=mean(m01)+(c*sd(m01))
  lim2=mean(m02)+(c*sd(m02))
  lim3=mean(m03)+(c*sd(m03))

#  if(diag.plot==T){

 #   if(type.plot=="all"){
  #    abline(h=lim1,lty=2,lwd=2)
   #   abline(h=lim2,lty=2,lwd=2)
    #  abline(h=lim3,lty=2,lwd=2)
    #}

    #if(type.plot=="rp"){
     # abline(h=lim1,lty=2,lwd=2)
    #}

    #if(type.plot=="smp"){
    #abline(h=lim2,lty=2,lwd=2)
    #}

    #if(type.plot=="evp"){
     #plot(m03,ylab="M(0)",pch=20)
     #abline(h=lim3,lty=2,lwd=2,pch=20)
    #}

  #}


  if(diag.plot==T){
    oldpar=par(mfrow=c(1,3))
    if(type.plot=="all"){
      on.exit(par(oldpar))
      plot1=plot(m01,ylab="M(0)",pch=20)
      abline(h=lim1,lty=2,lwd=2)
      on.exit(plot1)
      plot2= plot(m02,ylab="M(0)",pch=20)
      abline(h=lim2,lty=2,lwd=2)
      on.exit(plot2)
      plot3=plot(m03,ylab="M(0)",pch=20)
      abline(h=lim3,lty=2,lwd=2)
      on.exit(plot3)
      }

    if(type.plot=="rp"){
      plot1=plot(m01,ylab="M(0)",pch=20)
      on.exit(plot1)
    }

    if(type.plot=="smp"){
      plot2= plot(m02,ylab="M(0)",pch=20)
      abline(h=lim2,lty=2,lwd=2)
      on.exit(plot2)
    }

    if(type.plot=="evp"){
      plot3=plot(m03,ylab="M(0)",pch=20)
      abline(h=lim3,lty=2,lwd=2)
      on.exit(plot3)
    }

  }



  at1=at2=at3=NULL

  at1[m01>lim1]="atypical obs"
  at1[m01<=lim1]="normal obs"

  at2[m02>lim2]="atypical obs"
  at2[m02<=lim2]="normal obs"

  at3[m03>lim3]="atypical obs"
  at3[m03<=lim3]="normal obs"

  data.frame1=data.frame(at1,m01)
  data.frame2=data.frame(at2,m02)
  data.frame3=data.frame(at3,m03)


  s=list(Qwrp=Qw1,Qwsmp=Qw2,Qwevp=Qw3,respper=data.frame1,smper=data.frame2,expvper=data.frame3,limrp=lim1,limsmp=lim2,limevp=lim3)
  class(s)="LocInSCL"
  return(s)

}


