derivcormatrix=function(coords,phi,kappa=0,cov.model="exponential"){

  if(sum(is.na(coords)) > 0) stop("There are some NA values in the specified coordinates")
  if (!is.numeric(coords)& !is.data.frame(coords)) stop(" 2D coordinates must be a numeric matrix or data.frame")
  if (!is.matrix(coords)) coords=as.matrix(coords)
  if (ncol(as.matrix(coords)) !=2) stop("2D coordinates must be specified")

  if(!is.numeric(phi)) stop("phi must be a real number in [0,Inf)")
  if(!is.numeric(kappa)) stop("kappa must be a real number in [0,Inf)")
  if(phi<0)            stop("phi must be a real number in [0,Inf)")
  if(kappa<0) stop("kappa must be a real number in [0,Inf)")

  if(cov.model!="matern" & cov.model!="exponential" & cov.model!="gaussian" & cov.model!="spherical" &
     cov.model!="powered.exponential" & cov.model!="stable" & cov.model!="cauchy") {
    stop('Valid covariance structures are matern, exponential, gaussian, spherical,
         powered.exponential,stable, cauchy')
  }

  out=derivater(coords=coords,phi=phi,kappa=kappa,type=cov.model)

  s=list(H=out$H,d1=out$d1,d2=out$d2)

  return(s)
}
