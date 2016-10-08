\name{Seminaive}
\alias{Seminaive}

\title{Seminaive algorithm for spatial censored prediction.

}
\description{

This function executes the seminaive algorithm proposed by Schelin et al. (2014)}
\usage{
Seminaive(data, y.col, coords.col, covar, covar.col, copred,cov.model = "exponential",
thetaini, fix.nugget = T, nugget,kappa = 0, cons, MaxIter, cc, cutof, trend)
}

\arguments{
  \item{data}{

data.frame containing the coordinates, covariates and response variable.
}
  \item{y.col}{

(numeric) column of data.frame that corresponds to the response variable.
}
  \item{coords.col}{

(numeric) columns of data.frame that corresponds to the coordinates of the spatial data.
}
  \item{covar}{


(logical) indicates the presence of covariates in the spatial censored estimation (\code{FALSE}: without covariates, \code{TRUE}: with covariates).
}
  \item{covar.col}{

(numeric) columns of data.frame that corresponds to the covariates in the spatial censored linear model estimation.
}
  \item{copred}{
coordinates used in the prediction procedure.
}
  \item{cov.model}{
covariance model in the structure of covariance (see \code{cov.spatial} from \code{geoR}).
}
  \item{thetaini}{
initial values for the \eqn{\sigma^2} and \eqn{\phi} values in the covariance structure.
}
  \item{fix.nugget}{

(logical) it indicates if the \eqn{\tau^2} parameter must be fixed.
}
  \item{nugget}{
(numeric) values of the \eqn{\tau^2} parameter, if \code{fix.nugget=F}, this value corresponds to an initial value.
}
  \item{kappa}{

value of \eqn{\kappa} involved in some covariance functions.
}
  \item{cons}{
(vector) vector containing the \eqn{(c_1,c_2,c_3)} constants used in the convergence criterion for the algorithm (see Schedlin).
}
  \item{MaxIter}{
maximum of iterations for the algorithm.
}
  \item{cc}{
(binary vector) indicator of censure (1: censored, 0: observed)
}
  \item{cutof}{
(vector) limit of detection for censure ( rc: >cutoff, lc: <cutoff)
}
  \item{trend}{
it specifies the mean part of the model. See documentation of \code{trend.spatial} from \code{geoR} for further details. By default \code{"cte"}.
}
}
\details{
This function estimates and computes predictions following Schedlin et al. (2014). See reference.
}
\value{

\item{zk}{vector with observed and estimate censored observations by kriging prediction.}

\item{AIC}{AIC of the estimated model.}

\item{BIC}{BIC of the estimated model.}

\item{beta}{beta parameter for the mean structure.}

\item{theta}{vector of estimate parameters for the mean and covariance structure (\eqn{\beta,\sigma^2,\phi,\tau^2}).}

\item{predictions}{Predictions obtained for the seminaive algorithm.}

\item{sdpred}{Standard deviations of predictions.
}
\item{loglik}{log likelihood from the estimated model.}

}
\references{
Schelin, L. & Sjostedt-de Luna, S. (2014). Spatial prediction in the presence of left-censoring.
Computational Statistics and Data Analysis, 74.
}
\author{

Alejandro Ordonez <<ordonezjosealejandro@gmail.com>>, Victor H. Lachos <<hlachos@ime.unicamp.br>>
and Christian E. Galarza <<cgalarza88@gmail.com>>

Maintainer: Alejandro Ordonez <<ordonezjosealejandro@gmail.com>>
}

\seealso{
\code{\link{SAEMSCL}}
}

\examples{

\dontrun{

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
cutof=rep(obj$cutof,length(y[cc==1]))

###seminaive algorithm
r=Seminaive(data=data2,y.col=3,covar=T,coords.col=1:2,covar.col=4:5,cov.model="matern",
thetaini=c(.1,.2),fix.nugget=T,nugget=0,kappa=1.5,cons=c(0.1,2,0.5),MaxIter=100,
cc=obj$cc,cutof=cutof,copred=obj$coords1,trend=~V4+V5)

}
}


\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}
