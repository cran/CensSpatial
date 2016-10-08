rspacens=function(cov.pars,beta,x=as.matrix(rep(1,n)),coords,kappa=0,cens,n,n1,cov.model="exponential",cens.type){

  if (cens.type!='left'& cens.type!='right') stop('cens.type must be left or right')

  if(cov.model!="matern" & cov.model!="exponential" & cov.model!="gaussian" &cov.model!="spherical" &
     cov.model!="circular" & cov.model!="cubic" & cov.model!="wave" & cov.model!="linear" &
     cov.model!="power" &cov.model!="powered.exponential" &cov.model!="stable" & cov.model!="cauchy" &
     cov.model!="gencauchy" &cov.model!="gneiting" &cov.model!="gneiting.matern" &cov.model!="pure.nugget") {
    stop('cov.model should be one of matern, exponential, gaussian, spherical,
circular,cubic, wave, linear, power, powered.exponential, stable, cauchy, gencauchy,
         gneiting, gneiting.matern, pure.nugget')
  }

  if(!is.numeric(kappa)) stop("kappa must be a real number in [0,Inf)")
  if(kappa<0) stop("kappa must be a real number in [0,Inf)")

  if (!is.numeric(coords) & !is.data.frame(coords)) stop("2D coordinates must be a numeric matrix or data.frame")
  if (!is.matrix(coords)) coords=as.matrix(coords)
  if (!is.numeric(x)) stop("x must be a numeric matrix")
  if (!is.matrix(x)) x=as.matrix(x)
  if(!is.numeric(cens)) stop("cens must be a real number in [0,1]")
  if(cens > 1 | cens < 0) stop("cens must be a real number in [0,1]")
  if(n <= 1 |n%%1!=0) stop("n must be a positive integer value (greater than 1)")
  if(n1 <= 1 |n1%%1!=0) stop("n1 must be a positive integer value (greater than 1)")
  if(!is.numeric(beta)) stop("beta must be a numeric vector")
  if(length(beta)!=ncol(x)) stop("length of beta must be the same than the number of columns of x")
  if(!is.numeric(cov.pars)) stop("cov.pars must be a numeric vector")
  if(length(cov.pars)!=3) stop("the vector cov.pars=c(phi,sigma2,tau2) must be specified")

  ###validating NAS
  if(sum(is.na(coords)) > 0) stop("There are some NA values in coords")
  if(sum(is.na(x)) > 0) stop("There are some NA values in x")

  out=simcenspadata(cov.pars=cov.pars,beta=beta,x=x,coords=coords,kappa=kappa,cens=cens,n=n,n1=n1,type=cov.model,cens.type=cens.type)
  return(out)
}

