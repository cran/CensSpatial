\name{predSCL}
\alias{predSCL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Prediction for the SAEM algorithm for censored spatial data.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function uses the parameters estimates from SAEM to predict values at unknown locations through the MSE criterion assuming normal distribution.
}
\usage{
predSCL(xpred, coordspred, est)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xpred}{
%%     ~~Describe \code{xpred} here~~
values of the x design matrix for prediction coordinates.
}
  \item{coordspred}{
%%     ~~Describe \code{coordspred} here~~
points coordinates to be predicted.
}
  \item{est}{
%%     ~~Describe \code{est} here~~
object of the class SAEMSpatialCens (see \code{SAEMSCL} function).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function predicts using the Mean Square of error (MSE) criterion, that is, it takes the conditional
expectation \eqn{E(Y|X)} as the predictor that minimizes the MSE.
}
\value{
\item{prediction}{prediction value.}
\item{indpred}{indicator for the observed and predicted values (0:observed,1:predicted).}
\item{sdpred}{standard deviation for prediction.}
\item{coordspred}{points coordinates predicted.}
\item{coordsobs}{observed coordinates.}
}
\references{
DELYON, B., LAVIELLE, M.,ANDMOULI NES, E. (1999).  Convergence of a stochastic approximation version of the EM algorithm.Annals of Statistic-s27, 1, 94-128.

Diggle, P. & Ribeiro, P. (2007). Model-Based Geostatistics. Springer Series in Statistics.
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

###simulated coordinates
r1=sample(seq(1,30,length=400),n+n1)
r2=sample(seq(1,30,length=400),n+n1)

coords=cbind(r1,r2)### coordinates for estimation and prediction.

coords1=coords[1:n,]####coordinates used in estimation.

cov.ini=c(0.2,0.1)###initial values for phi and sigma2.

type="matern"
xtot<-cbind(1,runif((n+n1)),runif((n+n1),2,3))###X matrix for estimation and prediction.

xobs=xtot[1:n,]###X matrix for estimation.
beta=c(5,3,1)

###simulated data
obj=rspacens(cov.pars=c(3,.3,0),beta=beta,x=xtot,coords=coords,kappa=1.2,cens=0.25,
n=(n+n1),n1=n1,cov.model=type,cens.type="left")

data2=obj$datare
cc=obj$cc
y=obj$datare[,3]
coords=obj$datare[,1:2]
LI=rep(-Inf,length(y))
LS=rep(Inf,length(y))
LS[cc==1]=obj$cutof

#######SAEMSpatialCens object########

est=SAEMSCL(cc,y,cens.type="left",trend="cte",coords=coords,kappa=1.2,M=15,
perc=0.25,MaxIter=10,pc=0.2,cov.model="exponential",fix.nugget=TRUE,nugget=0,
inits.sigmae=cov.ini[2],inits.phi=cov.ini[1],search=TRUE,lower=0.00001,upper=50)


coordspred=obj$coords1
xpred=rep(1,dim(coordspred)[1])
h=predSCL(xpred,coordspred,est)
}

}

\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}
