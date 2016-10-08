localinfmeas=function(est,fix.nugget=T,diag.plot=T,type.plot="all",c=3){

  if(class(est)!="SAEMSpatialCens") stop("an object of the class SAEMSpatialCens must be provided")
  if(!is.logical(fix.nugget)) stop("fix.nugget must be TRUE or FALSE")
  if(!is.logical(diag.plot)) stop("diag.plot must be TRUE or FALSE")
  if(type.plot!="all" & type.plot!="rp" & type.plot!="smp" & type.plot!="evp"){
    stop("type.plot must be all,rp, smp,evp")
  }

  cov.model=est$type
  if(cov.model!="matern" & cov.model!="exponential" & cov.model!="gaussian" & cov.model!="spherical" &
     cov.model!="powered.exponential" & cov.model!="stable" & cov.model!="cauchy") {
    stop('Valid covariance structures are matern, exponential, gaussian, spherical,
         powered.exponential,stable, cauchy')
  }

  if(!is.numeric(c)) stop("the constant c must be a real number in [0,Inf)")
  if(c<0) stop("the constant c must be a real number in [0,Inf)")


  out=locinme(est=est,fix.nugget=fix.nugget,diag.plot=diag.plot,type.plot=type.plot,c=c)


  cat('\n\n')
  cat('---------------------------------------------------\n')
  cat('   Diagnostic in Spatial Censored Linear regression
         with Normal errors (SAEM estimation)   \n')
  cat('---------------------------------------------------\n')
  cat('\n')
  cat("Perturbation Schemes")
  cat('\n')
  cat('\n')
  cat("*Response perturbation \n")
  cat("*Scale matrix perturbation \n")
  cat("*Explanatory variable perturbation \n")
  cat('\n')

  return(out)

}
