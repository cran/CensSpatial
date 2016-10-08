Seminaive=function(data,y.col,coords.col,covar,covar.col,copred,cov.model="exponential",thetaini,fix.nugget=T,nugget,kappa=0,cons,MaxIter,cc,cutof,trend){

  if(!is.data.frame(data) & !is.numeric(data)) stop("data must be a data.frame or numeric matrix")
  if(!is.data.frame(data)) data=as.data.frame(data)
  if(!is.numeric(cutof)) stop ("cutof mus be a numeric vector")
  if(!is.numeric(thetaini)) stop ("thetaini must be a numeric vector thetaini=(sigma2,phi)")
  if(length(which(thetaini<0))>0) stop("all initials values for the covariance structure must be in [0,Inf)")
  if(length(thetaini)!=2) stop("thetaini must contains the initials values for sigma2 and phi")
  if(!is.numeric(kappa)) stop("kappa must be a real number in [0,Inf)")
  if(kappa<0) stop("kappa must be a real number in [0,Inf)")
  if(!is.numeric(nugget)) stop("kappa must be a real number in [0,Inf)")
  if(nugget<0) stop("nugget must be a real number in [0,Inf)")
  if(!is.logical(fix.nugget)) stop("fix.nugget must be TRUE or FALSE.")
  if(!is.logical(covar)) stop(" covar must be TRUE or FALSE.")
  if(nrow(data)!=nrow(cc)) stop ("data and cc does not have the same number of lines")
  if(length(which(cc==1))< length(cutof)) stop("Detections limits not coincide with the number of censored observations")
  if (trend!='cte'& trend!='1st' & trend!="2nd" & class(trend)!= "formula") stop('trend is not valid (see trend.spatial from geoR)')
  if(!is.numeric(y.col)) stop("a correct column for y (response) must be specified")
  if((y.col%%1)!=0) stop("a correct column for y (response) must be specified")
  if (sum(cc%in%c(0,1))< length(cc)) stop("The elements of the vector cc must be 0 or 1")
  if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
  if(sum(which(cons<0))>0) stop("all values in vector cons=(c1,c2,c3) must be positive")
  if(length(cons)!=3) stop("cons must contains (c1,c2,c3) constants for the algorithm")

  if(cov.model!="matern" & cov.model!="exponential" & cov.model!="gaussian" &cov.model!="spherical" &
     cov.model!="circular" & cov.model!="cubic" & cov.model!="wave" & cov.model!="linear" &
     cov.model!="power" &cov.model!="powered.exponential" &cov.model!="stable" & cov.model!="cauchy" &
     cov.model!="gencauchy" &cov.model!="gneiting" &cov.model!="gneiting.matern" &cov.model!="pure.nugget") {
    stop('cov.model should be one of matern, exponential, gaussian, spherical,
circular,cubic, wave, linear, power, powered.exponential, stable, cauchy, gencauchy,
gneiting, gneiting.matern, pure.nugget')
  }

  if(covar==T){
    if(!is.numeric(covar.col)) stop("a correct columns for  x (design matrix) must be specified")
    div1=covar.col%%1
    if(sum(which(div1!=0))>0) stop("a correct columns for  x (design matrix) must be specified")
  }

  ##validating NA
  if(sum(is.na(data)) > 0) stop("There are some NA values in data")
  if(sum(is.na(cc)) > 0) stop("There are some NA values in cc")
  if(sum(is.na(copred)) > 0) stop("There are some NA values in copred")
  if(sum(is.na(cutof)) > 0) stop("There are some NA values in cutof")

  ######Not data
  y = data[,y.col]
  coords = data[,coords.col]
  
  if((length(y) == 0) | (length(cc) == 0)|(length(coords) == 0)|(length(copred) == 0)|(length(cutof) == 0)){
    stop("All parameters must be provided")}

  if(covar==F){
    if(class(trend)=="formula") stop("trend must be cte,1stor 2nd when covar=F")
  }


  out=suppressWarnings(seminaive(data=data,y.col=y.col,covar=covar,coords.col=coords.col,covar.col=covar.col,cov.model=cov.model,thetaini=thetaini,fix.nugget=fix.nugget,nugget=nugget,kappa=kappa,cons=cons,MaxIter=MaxIter,cc=cc,cutof=cutof,copred=copred,trend=trend))



  if(trend=="cte"){
    trend1="Constant trend"
  }
  if(trend=="1st"){
    trend1="Linear function of its coordinates,

mu = beta0 + beta1*CoordX + beta2*CoordY"
  }
  if(trend=="2nd"){
    trend1="Linear function of its coordinates,

mu = beta0 + beta1*CoordX + beta2*CoordY + beta3*(CoordX)^2 +
   + beta4*(CoordY)^2 + beta5*(CoordX*CoordY)"
  }

  if(class(trend)=="formula"){
    trend1= "Linear trend,

mu = X*beta"
  }


  #Running the algorithm
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')
  out=out
  cat('\n\n')
  cat('---------------------------------------------------\n')
  cat('  Spatial Censored Linear regression with Normal errors (Seminaive estimation) \n')
  cat('---------------------------------------------------\n')
  cat('\n')
  cat("*Type of trend:",trend1)
  cat('\n')
  cat('\n')
  cat("*Covariance structure:",cov.model)
  cat('\n')
  cat('---------\n')
  cat('Estimates\n')
  cat('---------\n')
  cat('\n')
  trends= out$beta
  l = length(trends)

  lab = numeric(l+3)
  for (i in 1:l){ lab[i] = paste('beta ',i-1,sep='')}
  lab[l+1] = 'sigma2'
  lab[l+2] ='phi'
  lab[l+3] ='tau2'
  tab = round(cbind(out$theta),4)
  rownames(tab)=lab
  colnames(tab)=c("Seminaive est.")

  print(tab)
  cat('\n')
  cat('------------------------\n')
  cat('Model selection criteria\n')
  cat('------------------------\n')
  cat('\n')
  critFin <- c(out$loglik, out$AIC, out$BIC)
  critFin <- round(t(as.matrix(critFin)),digits=3)
  rownames(critFin) <- c("Value")
  colnames(critFin)=c("Loglik", "AIC", "BIC")
  print(critFin)
  cat('\n')


  return(out)

}