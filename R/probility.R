prob_opt = function(lower = rep(-Inf,ncol(sigma)),upper = rep(Inf,ncol(sigma)),mean = rep(0,ncol(sigma)),sigma,uselog2 = FALSE){
  p = ncol(sigma)
    #normal case
    if(p < 10){
      prob = ifelse(uselog2,
                    log2(pmvnorm(lower = lower,upper = upper,mean = mean,sigma = sigma)[1]),
                    pmvnorm(lower = lower,upper = upper,mean = mean,sigma = sigma)[1])
      if(prob < 0){
        prob = pmvn(lower = lower,upper = upper,mean = mean,sigma = sigma,uselog2 = uselog2)[[1]]
      }
    }else{
      prob = pmvn(lower = lower,upper = upper,mean = mean,sigma = sigma,uselog2 = uselog2)[[1]]
    }
  return(prob)
}
