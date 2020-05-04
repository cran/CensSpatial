\name{summary.seminaive}
\alias{summary.seminaive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of a seminaive object
%%  ~~function to do ... ~~
}
\description{
  \code{summary} method for class "seminaive".
}

\usage{
\method{summary}{seminaive}(object,\dots)
}

\arguments{
    \item{object}{
%%     ~~Describe \code{est} here~~
object of the class "seminaive" (see \code{\link{Seminaive}} function).
}
\item{\dots}{Additional arguments.}
   }

\value{
 \item{mean.str}{
Estimates for the mean structure parameters \eqn{\mathbf{beta}} for seminaive method.
}
 \item{var.str}{
Estimates for the variance structure parameters \eqn{\sigma^2, \phi} for seminaive method.
}


 \item{predictions}{
predictions for seminaive method.
}


}

\references{
Schelin, L. & Sjostedt-de Luna, S. (2014). Spatial prediction in the presence of left-censoring. Computational Statistics and Data Analysis, 74.
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
n<-10### sample size for estimation.
n1=1 ### number of observation used in the prediction.

###simulated coordinates.
r1=sample(seq(1,30,length=400),n+n1)
r2=sample(seq(1,30,length=400),n+n1)
coords=cbind(r1,r2)### total coordinates (used in estimation and prediction).
coords1=coords[1:n,]####coordinates used for estimation.

type="matern"### covariance structure.

xtot<-cbind(1,runif((n+n1)),runif((n+n1),2,3))## X matrix for estimation and prediction.
xobs=xtot[1:n,]## X matrix for estimation.

###simulated data.
obj=rspacens(cov.pars=c(3,.3,0),beta=c(5,3,1),x=xtot,coords=coords,kappa=1.2,
cens=0.25,n=(n+n1),n1=n1,cov.model=type,cens.type="left")

data2=obj$datare
data2[,4:5]=xobs[,-1]

cc=obj$cc
y=obj$datare[,3]
cutoff=rep(obj$cutoff,length(y[cc==1]))

###seminaive algorithm
r=Seminaive(data=data2,y.col=3,covar=TRUE,coords.col=1:2,covar.col=4:5,cov.model="matern",
thetaini=c(.1,.2),fix.nugget=TRUE,nugget=0,kappa=1.5,cons=c(0.1,2,0.5),MaxIter=100,
cc=obj$cc,cutoff=cutoff,copred=obj$coords1,trend=~V4+V5)

summary(r)

}

\donttest{
n<-200 ### sample size for estimation.
n1=100 ### number of observation used in the prediction.

###simulated coordinates.
r1=sample(seq(1,30,length=400),n+n1)
r2=sample(seq(1,30,length=400),n+n1)
coords=cbind(r1,r2)### total coordinates (used in estimation and prediction).
coords1=coords[1:n,]####coordinates used for estimation.

type="matern"### covariance structure.

xtot<-cbind(1,runif((n+n1)),runif((n+n1),2,3))## X matrix for estimation and prediction.
xobs=xtot[1:n,]## X matrix for estimation.

###simulated data.
obj=rspacens(cov.pars=c(3,.3,0),beta=c(5,3,1),x=xtot,coords=coords,kappa=1.2,
cens=0.25,n=(n+n1),n1=n1,cov.model=type,cens.type="left")

data2=obj$datare
data2[,4:5]=xobs[,-1]

cc=obj$cc
y=obj$datare[,3]
cutoff=rep(obj$cutoff,length(y[cc==1]))

###seminaive algorithm
r=Seminaive(data=data2,y.col=3,covar=TRUE,coords.col=1:2,covar.col=4:5,cov.model="matern",
thetaini=c(.1,.2),fix.nugget=TRUE,nugget=0,kappa=1.5,cons=c(0.1,2,0.5),MaxIter=100,
cc=obj$cc,cutoff=cutoff,copred=obj$coords1,trend=~V4+V5)

summary(r)

}


}

