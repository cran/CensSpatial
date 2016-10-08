#####SAEM_Spatial

SAEMSCL= function(cc, y, cens.type="left", trend="cte", LI=NULL, LS=NULL, x=NULL, coords, kappa=0, M = 20, perc = 0.25, MaxIter = 300, pc = 0.2,
                  cov.model="exponential", fix.nugget = TRUE, nugget, inits.sigmae, inits.phi,
                  search = F, lower, upper)
{

  m = length(y)
  if (trend!='cte'& trend!='1st' & trend!="2nd" & trend!="other") stop('trend must be cte, 1st,2nd or other')

  if(trend=="cte" | trend=="1st" | trend=="2nd"){

    if (!is.numeric(y)) stop("y must be a numeric vector")
    # if (!is.numeric(x)) stop("x must be a numeric matrix")
    #if (!is.numeric(x)) stop("x must be a numeric matrix")
    #if (!is.matrix(x)) x=as.matrix(x)
    if (!is.numeric(coords) & !is.data.frame(coords)) stop("2D coordinates must be a numeric matrix or data.frame")
    if (!is.matrix(coords)) coords=as.matrix(coords)
    #if (det(t(x)%*%x)==0) stop("the columns of x must be linearly independent")
    ## Verify error at parameters specification

    if (cens.type!='left'& cens.type!='right' & cens.type!="both") stop('cens.type must be left, right or both')


    if(cov.model!="matern" & cov.model!="exponential" & cov.model!="gaussian" &cov.model!="spherical" &
       cov.model!="circular" & cov.model!="cubic" & cov.model!="wave" & cov.model!="linear" &
       cov.model!="power" &cov.model!="powered.exponential" &cov.model!="stable" & cov.model!="cauchy" &
       cov.model!="gencauchy" &cov.model!="gneiting" &cov.model!="gneiting.matern" &cov.model!="pure.nugget") {
      stop('cov.model should be one of matern, exponential, gaussian, spherical,
circular,cubic, wave, linear, power, powered.exponential, stable, cauchy, gencauchy,
gneiting, gneiting.matern, pure.nugget')
    }


    #No data
    if( (length(y) == 0) | (length(cc) == 0)|(length(coords) == 0) ) stop("All parameters must be provided.")

    #Validating if exists NA's

    if (sum(cc%in%c(0,1))< length(cc)) stop("The elements of the vector cc must be 0 or 1")

    if(sum(is.na(y)) > 0) stop("There are some NA values in y")
    if(sum(is.na(coords)) > 0) stop("There are some NA values in coords")
    #if(sum(is.na(x)) > 0) stop("There are some NA values in x")
    if (sum(is.na(cc)) > 0) stop("There are some NA values in cc")


    #Validating dims data set
    if (ncol(as.matrix(y)) > 1) stop("y must have just one column")
    if (ncol(coords) !=2) stop("2D coordinates must be specified")
    if (ncol(as.matrix(cc)) > 1) stop("cc must have just one column")
    #if( length(y) != nrow(as.matrix(x)) ) stop("x does not have the same number of lines than y")
    if( length(y) != nrow(coords) ) stop("2D coordinates does not have the same number of lines than y")
    if( length(y) != length(cc) ) stop("cc does not have the same length than y")
    if(!is.numeric(MaxIter)) stop("MaxIter must be a positive integer value")
    if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
    if(!is.numeric(M)) stop("M must be a positive integer value")
    if(M <= 1 |M%%1!=0) stop("M must be a positive integer value (greater than 1)")
    if(!is.numeric(pc)) stop("pc must be a real number in [0,1]")
    if(pc > 1 | pc < 0) stop("pc must be a real number in [0,1]")
    if(!is.numeric(perc)) stop("perc must be a real number in [0,1)")
    if(perc >= 1 | perc < 0) stop("perc must be a real number in [0,1)")
    if(!is.numeric(inits.phi)) stop("Initial phi must be a real number in [0,Inf)")
    if(inits.phi<0)            stop("Initial phi must be a real number in [0,Inf)")
    if(!is.numeric(nugget)) stop("nugget must be a real number in [0,Inf)")
    if(nugget<0) stop("nugget must be a real number in [0,Inf)")
    if(!is.numeric(kappa)) stop("kappa must be a real number in [0,Inf)")
    if(kappa<0) stop("kappa must be a real number in [0,Inf)")
    if(!is.numeric(inits.sigmae)) stop("Initial sigma2 must be a real number in [0,Inf)")
    if(inits.sigmae<0) stop("Initial sigmae must be a real number in [0,Inf)")
    if(!is.logical(fix.nugget)) stop("fix.nugget must be TRUE or FALSE.")
    if(!is.logical(search)) stop("fix.nugget must be TRUE or FALSE.")

    if(cens.type=="left"|cens.type=="right"){
      if(!is.null(LI)|!is.null(LS)){ stop("The censure limits was specified by the indicator variable cc (not necessary to use this arguments when cens.type is left or right)")
      }

      out <-suppressWarnings(SAEM_Spatial(cc=cc,y=y,cens.type=cens.type,trend=trend,coords=coords,kappa=kappa,M=M,perc=perc,MaxIter=MaxIter,pc=pc,type.S=cov.model,
                                          fix.nugget=fix.nugget,nugget=nugget,inits.sigmae=inits.sigmae,inits.phi=inits.phi,search=search,lower=lower,upper=upper))

    }else{
      if(is.null(LI) | is.null(LS)) stop ("LI and LS must be specified when cens.type= both")
      if (!is.numeric(LI)) stop("LI must be a numeric vector")
      if (!is.numeric(LS)) stop("LS must be a numeric matrix")
      if( (length(LI) == 0) | (length(LS) == 0)) stop("LI and LS parameters must be provided in the presence of both types of censure.")
      if(sum(is.na(LS)) > 0) stop("There are some NA values in LI")
      if(sum(is.na(LI)) > 0) stop("There are some NA values in LS")
      if( length(y) != length(LI) ) stop("LI does not have the same number of lines than y")
      if( length(y) != length(LS) ) stop("LS does not have the same number of lines than y")

      out <-suppressWarnings(SAEM_Spatial(cc=cc,y=y,cens.type=cens.type,trend=trend,LI=LI,LS=LS,coords=coords,kappa=kappa,M=M,perc=perc,MaxIter=MaxIter,pc=pc,type.S=cov.model,
                                          fix.nugget=fix.nugget,nugget=nugget,inits.sigmae=inits.sigmae,inits.phi=inits.phi,search=search,lower=lower,upper=upper))
    }




  }



  if(trend=="other"){
    if(is.null(x)) stop("the trend matrix x must be specified")
    if (!is.numeric(y)) stop("y must be a numeric vector")
    if (!is.numeric(x)) stop("x must be a numeric matrix")
    if (!is.matrix(x)) x=as.matrix(x)
    if (!is.numeric(coords) & !is.data.frame(coords)) stop("2D coordinates must be a numeric matrix or data.frame")
    if (!is.matrix(coords)) coords=as.matrix(coords)
    if (det(t(x)%*%x)==0) stop("the columns of x must be linearly independent")
    ## Verify error at parameters specification

    if (cens.type!='left'& cens.type!='right' & cens.type!="both") stop('cens.type must be left, right or both')
    if (trend!='cte'& trend!='1st' & trend!="2nd" & trend!="other") stop('trend must be cte, 1st or 2nd or other')

    if(cov.model!="matern" & cov.model!="exponential" & cov.model!="gaussian" &cov.model!="spherical" &
       cov.model!="circular" & cov.model!="cubic" & cov.model!="wave" & cov.model!="linear" &
       cov.model!="power" &cov.model!="powered.exponential" &cov.model!="stable" & cov.model!="cauchy" &
       cov.model!="gencauchy" &cov.model!="gneiting" &cov.model!="gneiting.matern" &cov.model!="pure.nugget") {
      stop('cov.model should be one of matern, exponential, gaussian, spherical,
circular,cubic, wave, linear, power, powered.exponential, stable, cauchy, gencauchy,
           gneiting, gneiting.matern, pure.nugget')
    }
    #Validating LI and LS
    if(cens.type=="both"){
      if (!is.numeric(LI)) stop("LI must be a numeric vector")
      if (!is.numeric(LS)) stop("LS must be a numeric matrix")
      if( (length(LI) == 0) | (length(LS) == 0)) stop("LI and LS parameters must be provided in the presence of both types of censure.")
      if(sum(is.na(LS)) > 0) stop("There are some NA values in LI")
      if(sum(is.na(LI)) > 0) stop("There are some NA values in LS")
      if( length(y) != length(LI) ) stop("LI does not have the same number of lines than y")
      if( length(y) != length(LS) ) stop("LS does not have the same number of lines than y")
    }

    #No data
    if( (length(x) == 0) | (length(y) == 0) | (length(cc) == 0)|(length(coords) == 0) ) stop("All parameters must be provided.")

    #Validating if exists NA's

    if (sum(cc%in%c(0,1))< length(cc)) stop("The elements of the vector cc must be 0 or 1")

    if(sum(is.na(y)) > 0) stop("There are some NA values in y")
    if(sum(is.na(coords)) > 0) stop("There are some NA values in coords")
    if(sum(is.na(x)) > 0) stop("There are some NA values in x")
    if (sum(is.na(cc)) > 0) stop("There are some NA values in cc")


    #Validating dims data set
    if (ncol(as.matrix(y)) > 1) stop("y must have just one column")
    if (ncol(coords) !=2) stop("2D coordinates must be specified")
    if (ncol(as.matrix(cc)) > 1) stop("cc must have just one column")
    if( length(y) != nrow(as.matrix(x)) ) stop("x does not have the same number of lines than y")
    if( length(y) != length(cc) ) stop("cc does not have the same length than y")
    if( length(y) != nrow(coords) ) stop("2D coordinates does not have the same number of lines than y")
    if(!is.numeric(MaxIter)) stop("MaxIter must be a positive integer value")
    if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
    if(!is.numeric(M)) stop("M must be a positive integer value")
    if(M <= 1 |M%%1!=0) stop("M must be a positive integer value (greater than 1)")
    if(!is.numeric(pc)) stop("pc must be a real number in [0,1]")
    if(pc > 1 | pc < 0) stop("pc must be a real number in [0,1]")
    if(!is.numeric(perc)) stop("perc must be a real number in [0,1)")
    if(perc >= 1 | perc < 0) stop("perc must be a real number in [0,1)")
    if(!is.numeric(inits.phi)) stop("Initial phi must be a real number in [0,Inf)")
    if(inits.phi<0)            stop("Initial phi must be a real number in [0,Inf)")
    if(!is.numeric(nugget)) stop("nugget must be a real number in [0,Inf)")
    if(nugget<0) stop("nugget must be a real number in [0,Inf)")
    if(kappa<0) stop("kappa must be a real number in [0,Inf)")
    if(!is.numeric(inits.sigmae)) stop("Initial sigma2 must be a real number in [0,Inf)")
    if(inits.sigmae<0) stop("Initial sigmae must be a real number in [0,Inf)")
    if(!is.logical(fix.nugget)) stop("fix.nugget must be TRUE or FALSE.")
    if(!is.logical(search)) stop("search must be TRUE or FALSE.")


    if(cens.type=="left"|cens.type=="right"){
      if(!is.null(LI)|!is.null(LS)) stop("The censure limits was specified by the indicator variable cc (not necessary to use this arguments when cens.type is left or right)")
      out <-suppressWarnings(SAEM_Spatial(cc=cc,y=y,cens.type=cens.type,trend=trend,x=x,coords=coords,kappa=kappa,M=M,perc=perc,MaxIter=MaxIter,pc=pc,type.S=cov.model,
                                          fix.nugget=fix.nugget,nugget=nugget,inits.sigmae=inits.sigmae,inits.phi=inits.phi,search=search,lower=lower,upper=upper))
    }else{

      if (!is.numeric(LI)) stop("LI must be a numeric vector")
      if (!is.numeric(LS)) stop("LS must be a numeric matrix")
      if( (length(LI) == 0) | (length(LS) == 0)) stop("LI and LS parameters must be provided in the presence of both types of censure.")
      if(sum(is.na(LS)) > 0) stop("There are some NA values in LI")
      if(sum(is.na(LI)) > 0) stop("There are some NA values in LS")
      if( length(y) != length(LI) ) stop("LI does not have the same number of lines than y")
      if( length(y) != length(LS) ) stop("LS does not have the same number of lines than y")

      out <-suppressWarnings(SAEM_Spatial(cc=cc,y=y,cens.type=cens.type,trend=trend,LI=LI,LS=LS,x=x,coords=coords,kappa=kappa,M=M,perc=perc,MaxIter=MaxIter,pc=pc,type.S=cov.model,
                                          fix.nugget=fix.nugget,nugget=nugget,inits.sigmae=inits.sigmae,inits.phi=inits.phi,search=search,lower=lower,upper=upper))
    }

  }


  if(search==T){
    if(is.null(lower) | is.null(upper)) stop("lower and upper search limits must be specified")

    if(fix.nugget==T){
      if(length(lower)!=1 | length(upper)!=1) stop("specify a correct interval for phi in the real line")
    }

    if(fix.nugget==F){
      if(length(lower)!=2 | length(upper)!=2) stop("specify correct upper and lower limits for phi and tau2")
    }
  }

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

  if(trend=="other"){
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
  cat('  Spatial Censored Linear regression with Normal errors (SAEM estimation) \n')
  cat('---------------------------------------------------\n')
  cat('\n')
  cat("*Type of trend:",trend1)
  cat('\n')
  cat('\n')
  cat("*Covariance structure:",out$type)
  cat('\n')
  cat('---------\n')
  cat('Estimates\n')
  cat('---------\n')
  cat('\n')
  trends=out$X
  l = ncol(trends)
  if(fix.nugget){
    lab = numeric(l+2)
    for (i in 1:l) lab[i] = paste('beta ',i-1,sep='')
    lab[l+1] = 'sigma2'
    lab[l+2] ='phi'
    tab = round(cbind(out$theta,out$ep),4)
    rownames(tab)=t(lab)
    colnames(tab)="Estimated"
  }else{
    lab = numeric(l+3)
    for (i in 1:l) lab[i] = paste('beta ',i-1,sep='')
    lab[l+1] = 'sigma2'
    lab[l+2] ='phi'
    lab[l+3] ='tau2'
    tab = round(cbind(out$theta,out$ep),4)
    rownames(tab)=lab
    colnames(tab)="Estimated"
  }
  print(tab)
  cat('\n')
  cat('------------------------\n')
  cat('Model selection criteria\n')
  cat('------------------------\n')
  cat('\n')
  critFin <- c(out$loglik, out$AIC, out$BIC, out$AICcorr)
  critFin <- round(t(as.matrix(critFin)),digits=3)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","AICcorr"))
  print(critFin)
  cat('\n')
  cat('-------\n')
  cat('Details\n')
  cat('-------\n')
  cat('\n')
  cat('Type of censoring =',cens.type)
  cat('\n')
  if (sum(cc)>0) {
    cat("Convergence reached? =",(out$iter < MaxIter))
    cat('\n')
    cat('Iterations =',out$iter,"/",MaxIter)
    cat('\n')
    cat('MC sample =',M)
    cat('\n')
    cat('Cut point =',pc)
    cat('\n')
  }






  obj.out <- list(beta = out$beta1, sigma2 = out$sigmae, phi = out$phi, nugget = out$tau2, Theta=out$Theta, loglik=out$loglik,
                  AIC=out$AIC, BIC=out$BIC, AICcorr=out$AICcorr,X=out$X, Psi=out$Psi,trend=out$trend,
                  theta = out$theta, uyy = out$yy,uy=out$uy,cc=out$cc,type=out$type,kappa=out$kappa,coords=out$coords,iterations=out$iter,timex=out$timex,fitted=out$fitted)

  class(obj.out) <- "SAEMSpatialCens"

  return(invisible(obj.out))


}



#est2= SAEMSCL(cc,y,cens.type="left",trend="other",x=xobs,coords=coords,kappa=0.3,M=15,perc=0.25,MaxIter=4,pc=0.2,cov.model="spherical",
#                   fix.nugget=T,nugget=0,inits.sigmae=cov.ini[2],inits.phi=cov.ini[1],search=T,lower=0.00001,upper=50)








