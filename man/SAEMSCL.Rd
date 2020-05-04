\name{SAEMSCL}
\alias{SAEMSCL}
\title{SAEM Algorithm estimation for censored spatial data.

}
\description{
It estimates the parameters for a linear spatial model with censored observations
}


\usage{
SAEMSCL(cc, y, cens.type="left", trend = "cte", LI = NULL, LS = NULL, x = NULL, coords,
kappa = 0, M = 20, perc = 0.25, MaxIter = 300, pc = 0.2, cov.model = "exponential",
fix.nugget = TRUE, nugget, inits.sigmae, inits.phi, search = FALSE, lower, upper)
}

\arguments{
  \item{cc}{
(binary vector) indicator of censure (1: censored observation 0: observed).
}
  \item{y}{

(vector) corresponds to  response variable.
}
  \item{cens.type}{
type of censure (\code{"left"}:left or \code{"right"}:right).
}
  \item{trend}{
linear trends options: \code{"cte"}, \code{"1st"}, \code{"2nd"} and \code{"other"}, the three first are defined
like in \code{geoR}, if trend=\code{"other"}, x (design matrix) must be defined.
}
  \item{LI}{
(vector) lower limit, if cens.type=\code{"both"}, LI must be provided, if cens.type=\code{"left"}
or \code{"right"} LI and LS are defined by the function through the indicator of censure cc.
}
  \item{LS}{
(vector) upper limit, if \code{cens.type="both"}, LS must be provided, if \code{cens.type="left"}
or "right" LI and LS are defined by the function through the indicator of censure cc.

}
  \item{x}{
design matrix.
}
  \item{coords}{
corresponds to the coordinates of the spatial data (2D coordinates).
}
  \item{kappa}{
value of kappa used in some covariance functions.
}
  \item{M}{
number of montecarlo samples for stochastic aproximation.
}
  \item{perc}{
percentage of burn-in on the Monte Carlo sample. Default=0.25.
}
  \item{MaxIter}{
maximum of iterations for the algorithm.
}
  \item{pc}{
percentage of initial iterations of the SAEM algorithm.
(Default=0.2).
}
  \item{cov.model}{
covariance Structure (see, \code{cov.spatial} from \code{geoR}).
}
  \item{fix.nugget}{
(logical) indicates if the \eqn{\tau^2} parameter must be fixed.
}
  \item{nugget}{
if \code{fix.nugget=TRUE}, the algorithm just estimates \eqn{\beta}, \eqn{\sigma^2}, and \eqn{\phi}, and fixed \eqn{\tau^2} like nugget, else, \eqn{\tau^2} is estimated and nugget corresponds to initial value for \eqn{\tau^2}.

}
  \item{inits.sigmae}{

 corresponds to initial value for \eqn{\sigma^2}.
}
  \item{inits.phi}{

corresponds to initial value for \eqn{\phi} parameter.
}
  \item{search}{

(logical) this argument gives bounds where the optim routine can find the solution that maximizes the Maximum likelihood expectation. If \code{search=F}, the optim routine will try to search the solutions for maximization in all the domain for \eqn{\phi} and \eqn{\tau^2} (if \code{fix.nugget=FALSE}). If \code{search=TRUE}, the optim routine search the solutions in a specific neighborhood. We recommended to use \code{search=F} (see details).

}
  \item{lower}{
(vector or numeric) lower bound from the optim solution. If \code{fix.nugget=T}, lower is numerical and corresponds to the lower bound for search the solution of the \eqn{\phi} parameter, if \code{fix.nugget=FALSE} lower is a vector and corresponds to the lower bounds for search the solution of \eqn{\phi} and \eqn{\tau2} that maximizes the Maximum Likelihood Expectation
(see details).
}
  \item{upper}{
(vector or numeric) upper bound from the optim solution. If \code{fix.nugget==T}, lower is numerical and corresponds to the lower bound for searching the solution of the phi parameter, if \code{fix.nugget==F}, lower is a vector and corresponds to the lower bounds for searching the solution for \eqn{\phi} and \eqn{\tau^2} parameters that maximizes the Maximum Likelihood Expectation}
}
\details{
The estimation process was computed via SAEM algorithm initially proposed by Deylon et. al.(1999). This is a stochastic approximation of the EM procedure. This procedure circunvent the heavy computational time involved in the MCEM procedure necessary for estimating phi and tau2 parameters
(when tau2 is not fixed) since there is not an analytical solution. The search interval was proposed because sometimes the maximization procedure used by optim function does not work for large intervals.
}
\value{
\item{beta}{estimated \eqn{\beta}.}
\item{sigma2}{estimated \eqn{\sigma^2}.}
\item{phi}{estimated \eqn{\phi}.}
\item{nugget}{estimated or fixed \eqn{\tau^2}.}
\item{Theta}{estimated parameters in all iterations (\eqn{\beta, \sigma^2, \phi}) or (\eqn{\beta, \sigma^2, \phi, \tau^2}) if \code{fix.nugget=F}.}
\item{loglik}{log likelihood for SAEM method.}
\item{AIC}{Akaike information criteria.}
\item{BIC}{Bayesian information criteria.}
\item{AICcorr}{corrected AIC by the number of parameters.}
\item{X}{design matrix.}
\item{Psi}{estimated covariance matrix.}
\item{theta}{final estimation of \eqn{\theta=(\beta, \sigma^2, \phi)} or \eqn{\theta=(\beta, \sigma^2, \phi, \tau^2)} if \code{fix.nugget=F}.}
\item{uy}{stochastic approximation of the first moment for the truncated normal distribution.}
\item{uyy}{stochastic approximation of the second moment for the truncated normal distribution.}
\item{cc}{indicator of censure (0:observed, 1: censored).}
\item{type}{covariance structure considered in the model.}
\item{kappa}{\eqn{\kappa} parameter for some covariance structures.}
\item{coords}{coordinates of the observed data.}
\item{iterations}{number of iterations needed to convergence.}
\item{fitted}{fitted values for the SAEM algortihm.}

}
\references{
DELYON, B., LAVIELLE, M.,ANDMOULI NES, E. (1999).  Convergence ofa stochastic approximation version of the EM algorithm.Annals of Statistic-s27, 1, 94-128.

Diggle, P. & Ribeiro, P. (2007). Model-Based Geostatistics. Springer Series in Statistics.

}
\author{
Alejandro Ordonez <<ordonezjosealejandro@gmail.com>>, Victor H. Lachos <<hlachos@ime.unicamp.br>>
and Christian E. Galarza <<cgalarza88@gmail.com>>

Maintainer: Alejandro Ordonez <<ordonezjosealejandro@gmail.com>>
}

\seealso{
\code{\link{localinfmeas}}, \code{\link{derivQfun}}
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
kappa=1.2,M=15,perc=0.25,MaxIter=5,pc=0.2,cov.model=type,
fix.nugget=TRUE,nugget=0,inits.sigmae=cov.ini[1],
inits.phi=cov.ini[2],search=TRUE,lower=0.00001,upper=100)

summary(est)
}


\donttest{
n<-200 ### sample size for estimation.
n1=100 ### number of observation used in the prediction.

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
}
\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}
