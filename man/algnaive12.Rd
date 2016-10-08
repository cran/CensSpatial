\name{algnaive12}
\alias{algnaive12}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Naive 1 and Naive 2 method for spatial prediction.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function performs spatial censored estimation and prediction for left and right censure through the Naive 1 and Naive 2 methods.
}
\usage{
algnaive12(data, cc, copred, thetaini, y.col = 3,coords.col = 1:2,covar=F, covar.col,
fix.nugget = T, nugget, kappa = 0, cutoff, cov.model = "exponential", trend)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data}{
%%     ~~Describe \code{data} here~~
data.frame containing the coordinates, covariates and the response variable (in any order).
}

\item{cc}{
%%     ~~Describe \code{cc} here~~
(binary vector) indicator of censure (1: censored observation 0: observed).
}

 \item{copred}{
%%     ~~Describe \code{copred} here~~
coordinates used in the prediction procedure.
}

 \item{thetaini}{
%%     ~~Describe \code{thetaini} here~~
initial values for the \eqn{\sigma^2} and \eqn{\phi} values in the covariance structure.
}

  \item{y.col}{
%%     ~~Describe \code{y.col} here~~
(numeric) column of data.frame that corresponds to the response variable.
}
  \item{coords.col}{
%%     ~~Describe \code{coords.col} here~~
(numeric) columns of data.frame that corresponds to the coordinates of the spatial data.
}
  \item{covar}{
%%     ~~Describe \code{covar} here~~

(logical) indicates the presence of covariates in the spatial censored estimation
(\code{FALSE} :without covariates, \code{TRUE} :with covariates).
}
  \item{covar.col}{
%%     ~~Describe \code{covar.col} here~~
(numeric) columns of data.frame that corresponds to the covariates in the spatial censored linear model estimation.
}


  \item{fix.nugget}{
%%     ~~Describe \code{fix.nugget} here~~
(logical) it indicates if the \eqn{\tau^2} parameter must be fixed.
}
  \item{nugget}{
%%     ~~Describe \code{nugget} here~~
(numeric) values of the \eqn{\tau^2} parameter, if \code{fix.nugget=F} this value corresponds to an initial value.
}
  \item{kappa}{
%%     ~~Describe \code{kappa} here~~
value of \eqn{\kappa} used in some covariance functions.
}

  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
(vector) Limit of censure detection ( rc:>cutoff, lc:<cutoff).
}

  \item{cov.model}{
%%     ~~Describe \code{cov.model} here~~
structure of covariance (see \code{cov.spatial} from \code{geoR}).
}

  \item{trend}{
%%     ~~Describe \code{trend} here~~
it specifies the mean part of the model. See documentation of trend.spatial from geoR for further details. By default it takes \code{"cte"}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The Naive 1 and Naive 2 are computed as in Schelin (2014). The naive 1 replaces the censored observations by the limit of detection (LD) and it performs estimation and prediction with this data. Instead of 1, the naive 2 replaces the censored observations by LD/2.

}
\value{

\item{beta1}{beta parameter for the mean structure in the Naive 1 method.}

\item{beta2}{beta parameter for the mean structure in the Naive 2 method.}

\item{theta1}{vector of estimate parameter for the mean and covariance structure (\eqn{\beta, \sigma^2, \phi, \tau^2}) in the
Naive 1 method.}

\item{theta2}{vector of estimate parameter for the mean and covariance structure (\eqn{\beta, \sigma^2, \phi, \tau^2}) in the
Naive 2 method.}

\item{predictions1}{predictions obtained for the Naive 1 method.}

\item{predictions2}{predictions obtained for the Naive 2 method.}

\item{AIC1}{AIC of the estimated model in the Naive 1 method.}

\item{AIC2}{AIC of the estimated model in the Naive 2 method.}

\item{BIC1}{BIC of the estimated model in the Naive 1 method.}

\item{BIC2}{BIC of the estimated model in the Naive 2 method.}

\item{loglik1}{log likelihood for the estimated model in the Naive 1 method.}

\item{loglik2}{log likelihood for the estimated model in the Naive 2 method.}

\item{sdpred1}{standard deviations of predictions in the Naive 1 method.
}

\item{sdpred2}{standard deviations of predictions in the Naive 2 method.
}


}
\references{
Schelin, L. & Sjostedt-de Luna, S. (2014). Spatial prediction in the presence of left-censoring.
Computational Statistics and Data Analysis, 74.
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

\dontrun{

n<-200 ### sample size for estimation.
n1=100 ### number of observation used for prediction.

###simulated coordinates
n<-200 ### sample size for estimation.
n1=100 ### number of observation used in the prediction.

###simulated coordinates
r1=sample(seq(1,30,length=400),n+n1)
r2=sample(seq(1,30,length=400),n+n1)
coords=cbind(r1,r2)### total coordinates (used in estimation and prediction).

coords1=coords[1:n,]####coordinates used for estimation.

type="matern"### covariance structure.

xtot<-cbind(1,runif((n+n1)),runif((n+n1),2,3))## X matrix for estimation and prediction.
xobs=xtot[1:n,]## X matrix for estimation.

###simulated data
obj=rspacens(cov.pars=c(3,.3,0),beta=c(5,3,1),x=xtot,coords=coords,kappa=1.2,
cens=0.25,n=(n+n1),n1=n1,cov.model=type,cens.type="left")

data2=obj$datare
data2[,4:5]=xobs[,-1]

cc=obj$cc
y=obj$datare[,3]
cutoff=rep(obj$cutoff,length(y[cc==1]))

aux2=algnaive12(data=data2,cc=obj$cc,covar=T,covar.col=4:5,copred=obj$coords1,
thetaini=c(.1,..2),y.col=3,coords.col=1:2,fix.nugget=T,nugget=0,kappa=1.2,
cutoff=cutoff,trend=~V4+V5,cov.model=type)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}
