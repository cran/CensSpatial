\name{summary.SAEMSpatialCens}
\alias{summary.SAEMSpatialCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of a SAEMSpatialCens object.
%%  ~~function to do ... ~~
}
\description{
  \code{summary} method for class "SAEMSpatialCens".
}

\usage{
\method{summary}{SAEMSpatialCens}(object,\dots)
}

\arguments{
    \item{object}{
%%     ~~Describe \code{est} here~~
object of the class "SAEMSpatialCens" (see \code{\link{SAEMSCL}} function).
}

\item{\dots}{Additional arguments.}
   }

\value{
 \item{mean.str}{
Estimates for the mean structure parameters \eqn{\mathbf{beta}} for SAEMSCL method.
}
 \item{var.str}{
Estimates for the variance structure parameters \eqn{\sigma^2, \phi} for SAEMSCL method.
}


}

\references{
DELYON, B., LAVIELLE, M.,ANDMOULI NES, E. (1999). Convergence of a stochastic approximation version of the EM algorithm.Annals of Statistic-s27, 1, 94-128.

Diggle, P. & Ribeiro, P. (2007). Model-Based Geostatistics. Springer Series in Statistics.
}

\author{
%%  ~~who you are~~
Alejandro Ordonez <<ordonezjosealejandro@gmail.com>>, Victor H. Lachos <<hlachos@ime.unicamp.br>>
and Christian E. Galarza <<cgalarza88@gmail.com>>

Maintainer: Alejandro Ordonez <<ordonezjosealejandro@gmail.com>>
}


\seealso{
\code{\link{SAEMSCL}}
}

\examples{

\dontshow{
n<-10 ### sample size for estimation.
n1=1 ### number of observation used in the prediction.

###simulated coordinates
r1=sample(seq(1,30,length=400),n+n1)
r2=sample(seq(1,30,length=400),n+n1)
coords=cbind(r1,r2)

coords1=coords[1:n,]

type="matern"
#xtot<-cbind(1,runif((n+n1)),runif((n+n1),2,3))
xtot=as.matrix(rep(1,(n+n1)))
xobs=xtot[1:n,]
beta=5
#beta=c(5,3,1)

###simulated data
obj=rspacens(cov.pars=c(3,.3,0),beta=beta,x=xtot,coords=coords,kappa=1.2,
cens=0.25,n=(n+n1),n1=n1,cov.model=type,cens.type="left")

data2=obj$datare
cc=obj$cc
y=obj$datare[,3]
coords=obj$datare[,1:2]
##initials values obtained from variofit.
cov.ini=c(0.13,0.86)

est=SAEMSCL(cc,y,cens.type="left",trend="cte",coords=coords,
kappa=1.2,M=15,perc=0.25,MaxIter=10,pc=0.2,cov.model=type,
fix.nugget=TRUE,nugget=0,inits.sigmae=cov.ini[1],
inits.phi=cov.ini[2],search=TRUE,lower=0.00001,upper=100)

summary(est)


}

\donttest{
n<-200 ### sample size for estimation.
n1=50 ### number of observation used in the prediction.

###simulated coordinates
r1=sample(seq(1,30,length=400),n+n1)
r2=sample(seq(1,30,length=400),n+n1)
coords=cbind(r1,r2)

coords1=coords[1:n,]

type="matern"
#xtot<-cbind(1,runif((n+n1)),runif((n+n1),2,3))
xtot=as.matrix(rep(1,(n+n1)))
xobs=xtot[1:n,]
beta=5
#beta=c(5,3,1)

###simulated data
obj=rspacens(cov.pars=c(3,.3,0),beta=beta,x=xtot,coords=coords,kappa=1.2,
cens=0.25,n=(n+n1),n1=n1,cov.model=type,cens.type="left")

data2=obj$datare
cc=obj$cc
y=obj$datare[,3]
coords=obj$datare[,1:2]

est=SAEMSCL(cc,y,cens.type="left",trend="cte",coords=coords,
kappa=1.2,M=15,perc=0.25,MaxIter=10,pc=0.2,cov.model=type,
fix.nugget=TRUE,nugget=0,inits.sigmae=cov.ini[1],
inits.phi=cov.ini[2],search=TRUE,lower=0.00001,upper=100)

summary(est)

}


}

