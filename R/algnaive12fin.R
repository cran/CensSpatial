algnaive12=function(data,cc,copred,thetaini,y.col=3,coords.col=1:2,covar=F,covar.col,fix.nugget=T,nugget,kappa=0,cutoff,cov.model="exponential",trend){

  if(!is.data.frame(data) & !is.numeric(data)) stop("data must be a data.frame or numeric matrix")
  if(!is.data.frame(data)) data=as.data.frame(data)
  if(!is.numeric(cutoff)) stop ("cutoff mus be a numeric vector")
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
  if(length(which(cc==1))< length(cutoff)) stop("Detections limits not coincide with the number of censored observations")
  if (trend!='cte'& trend!='1st' & trend!="2nd" & class(trend)!= "formula") stop('trend is not valid (see trend.spatial from geoR)')
  if(!is.numeric(y.col)) stop("a correct column for y (response) must be specified")
  if((y.col%%1)!=0) stop("a correct column for y (response) must be specified")
  if (sum(cc%in%c(0,1))< length(cc)) stop("The elements of the vector cc must be 0 or 1")

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
  if(sum(is.na(cutoff)) > 0) stop("There are some NA values in cutoff")

  ######Not data
  y = data[,y.col]
  coords = data[,coords.col]
                
  if((length(y) == 0) | (length(cc) == 0)|(length(coords) == 0)|(length(copred) == 0)|(length(cutoff) == 0)){
    stop("All parameters must be provided")}

  if(covar==F){
    if(class(trend)=="formula") stop("trend must be cte,1stor 2nd when covar=F")
  }


  out=naive12(data=data,cc=cc,covar=covar,covar.col=covar.col,copred=copred,thetaini=thetaini,y.col=y.col,coords.col=coords.col,fix.nugget=fix.nugget,nugget=nugget,kappa=kappa,cutof=cutoff,trend=trend)


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

  #out=naive12(data=data,cc=cc,covar=covar,covar.col=covar.col,copred=copred,thetaini=thetaini,y.col=y.col,coords.col=coords.col,fix.nugget=fix.nugget,nugget=nugget,kappa=kappa,cutof=cutoff,trend=trend)

  #Running the algorithm
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')
  out=out
  cat('\n\n')
  cat('---------------------------------------------------\n')
  cat('  Spatial Censored Linear regression with Normal errors (Naive 1 and Naive 2 estimation) \n')
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
  trends= out$beta1
  l = length(trends)

  lab = numeric(l+3)
  for (i in 1:l){ lab[i] = paste('beta ',i-1,sep='')}
  lab[l+1] = 'sigma2'
  lab[l+2] ='phi'
  lab[l+3] ='tau2'
  tab = round(cbind(out$theta1,out$theta2),4)
  rownames(tab)=lab
  colnames(tab)=c("Naive 1","Naive 2")

  print(tab)
  cat('\n')
  cat('------------------------\n')
  cat('Model selection criteria\n')
  cat('------------------------\n')
  cat('\n')
  critFin1 <- c(out$loglik1, out$AIC1, out$BIC1)
  critFin2 <- c(out$loglik2, out$AIC2, out$BIC2)
  critFin=rbind(critFin1,critFin2)
  critFin <- round(as.matrix(critFin),digits=3)
  rownames(critFin) <- c("Naive 1", "Naive 2")
  colnames(critFin)=c("Loglik", "AIC", "BIC")
  print(critFin)
  cat('\n')


  return(out)

}
