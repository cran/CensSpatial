summary.naive=function(object, ...){
  if(class(object)!="naive") stop ("An object of class naive must be provided")
  #Running the algorithm
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')
  object=object
  cat('\n\n')
  cat('---------------------------------------------------\n')
  cat('  Spatial Censored Linear regression with Normal errors (Naive 1 and Naive 2 estimation) \n')
  cat('---------------------------------------------------\n')
  cat('\n')
  cat("*Type of trend:",object$trend1)
  cat('\n')
  cat('\n')
  cat("*Covariance structure:",object$type)
  cat('\n')
  cat('---------\n')
  cat('Estimates\n')
  cat('---------\n')
  cat('\n')
  trends= object$beta1
  l = length(trends)

  lab = numeric(l+3)
  for (i in 1:l){ lab[i] = paste('beta ',i-1,sep='')}
  lab[l+1] = 'sigma2'
  lab[l+2] ='phi'
  lab[l+3] ='tau2'
  tab = round(cbind(object$theta1,object$theta2),4)
  rownames(tab)=lab
  colnames(tab)=c("Naive 1","Naive 2")

  print(tab)
  cat('\n')
  cat('------------------------\n')
  cat('Model selection criteria\n')
  cat('------------------------\n')
  cat('\n')
  critFin1 <- c(object$loglik1, object$AIC1, object$BIC1)
  critFin2 <- c(object$loglik2, object$AIC2, object$BIC2)
  critFin=rbind(critFin1,critFin2)
  critFin <- round(as.matrix(critFin),digits=3)
  rownames(critFin) <- c("Naive 1", "Naive 2")
  colnames(critFin)=c("Loglik", "AIC", "BIC")
  print(critFin)
  cat('\n')
  invisible(list(mean.str1=object$theta1[1:l],var.str1=object$theta1[(l+1):(l+2)],mean.str2=object$theta2[1:l],var.str2=object$theta2[(l+1):(l+2)],pred1=object$predictions1,pred2=object$predictions2))
  }

summary.seminaive=function(object, ...){
  if(class(object)!="seminaive") stop ("An object of class seminaive must be provided")
  #Running the algorithm
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')
  object=object
  cat('\n\n')
  cat('---------------------------------------------------\n')
  cat('  Spatial Censored Linear regression with Normal errors (Seminaive estimation) \n')
  cat('---------------------------------------------------\n')
  cat('\n')
  cat("*Type of trend:",object$trend1)
  cat('\n')
  cat('\n')
  cat("*Covariance structure:",object$type)
  cat('\n')
  cat('---------\n')
  cat('Estimates\n')
  cat('---------\n')
  cat('\n')
  trends= object$beta
  l = length(trends)

  lab = numeric(l+3)
  for (i in 1:l){ lab[i] = paste('beta ',i-1,sep='')}
  lab[l+1] = 'sigma2'
  lab[l+2] ='phi'
  lab[l+3] ='tau2'
  tab = round(cbind(object$theta),4)
  rownames(tab)=lab
  colnames(tab)=c("Seminaive est.")

  print(tab)
  cat('\n')
  cat('------------------------\n')
  cat('Model selection criteria\n')
  cat('------------------------\n')
  cat('\n')
  critFin <- c(object$loglik, object$AIC, object$BIC)
  critFin <- round(t(as.matrix(critFin)),digits=3)
  rownames(critFin) <- c("Value")
  colnames(critFin)=c("Loglik", "AIC", "BIC")
  print(critFin)
  cat('\n')

  invisible(list(mean.str=object$theta1[1:l],var.str=object$theta1[(l+1):(l+2)],pred=object$predictions))
  }


summary.SAEMSpatialCens=function(object, ...){
  if(class(object)!="SAEMSpatialCens") stop("An object of class SAEMSpatialCens must be provided")
  #Running the algorithm
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')

  object=object

  cat('\n\n')
  cat('---------------------------------------------------\n')
  cat('  Spatial Censored Linear regression with Normal errors (SAEM estimation) \n')
  cat('---------------------------------------------------\n')
  cat('\n')
  cat("*Type of trend:",object$trend1)
  cat('\n')
  cat('\n')
  cat("*Covariance structure:",object$type)
  cat('\n')
  cat('---------\n')
  cat('Estimates\n')
  cat('---------\n')
  cat('\n')
  trends=object$X
  l = ncol(trends)
  if(object$fix.nugget){
    lab = numeric(l+2)
    for (i in 1:l) lab[i] = paste('beta ',i-1,sep='')
    lab[l+1] = 'sigma2'
    lab[l+2] ='phi'
    tab = round(cbind(object$theta),4)
    rownames(tab)=t(lab)
    colnames(tab)="Estimated"
  }else{
    lab = numeric(l+3)
    for (i in 1:l) lab[i] = paste('beta ',i-1,sep='')
    lab[l+1] = 'sigma2'
    lab[l+2] ='phi'
    lab[l+3] ='tau2'
    tab = round(cbind(object$theta),4)
    rownames(tab)=lab
    colnames(tab)="Estimated"
  }
  print(tab)
  cat('\n')
  cat('------------------------\n')
  cat('Model selection criteria\n')
  cat('------------------------\n')
  cat('\n')
  critFin <- c(object$loglik, object$AIC, object$BIC, object$AICcorr)
  critFin <- round(t(as.matrix(critFin)),digits=3)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","AICcorr"))
  print(critFin)
  cat('\n')
  cat('-------\n')
  cat('Details\n')
  cat('-------\n')
  cat('\n')
  cat('Type of censoring =',object$cens.type)
  cat('\n')
  if (sum(object$cc)>0) {
    cat("Convergence reached? =",(object$iterations < object$MaxIter))
    cat('\n')
    cat('Iterations =',object$iterations,"/",object$MaxIter)
    cat('\n')
    cat('MC sample =',object$M)
    cat('\n')
    cat('Cut point =',object$pc)
    cat('\n')
  }

  invisible(list(mean.str=object$theta1[1:l],var.str=object$theta1[(l+1):(l+2)]))
  }

