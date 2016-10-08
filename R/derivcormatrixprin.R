
#########5.First and second derivates of some correlation matrix#########


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


  if(type=="powered.exponential"|type=="stable"){
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

