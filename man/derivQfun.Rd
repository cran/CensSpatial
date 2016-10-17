\name{derivQfun}
\alias{derivQfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum Likelihood Expectation (\eqn{logQ} function and its derivates)
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
It computes the \eqn{logQ} function, its derivates of first and second order and the inverse of the hessian matrix for the SAEM estimated parameters.
}
\usage{
derivQfun(est, fix.nugget = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est}{
%%     ~~Describe \code{est} here~~
object of the class "SAEMSpatialCens". See \code{SAEMSCL} function.
}
  \item{fix.nugget}{
%%     ~~Describe \code{fix.nugget} here~~
(logical) it indicates if the \eqn{\tau^2} parameter must be fixed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The \eqn{logQ} function refers to the logarithm of the Maximum likelihood conditional expectation, the first and second moments of the truncated normal distribution of censored data are involved in its computation.
}
\value{
\item{Qlogvalue}{value of the \eqn{logQ} function evaluated in the SAEM estimates.}

\item{gradQ}{gradient for the \eqn{logQ} function evaluated in the SAEM estimates.}

\item{HQ}{hessian Matrix for the \eqn{logQ} function evaluated in the SAEM estimates.}

\item{Qinv}{inverse of the negative Hessian matrix for the \eqn{logQ} function evaluated in the SAEM estimates.}
}
\references{
Diggle, P. & Ribeiro, P. (2007). Model-Based Geostatistics. Springer Series in Statistics.

Gradshtejn, I. S. & Ryzhik, I. M. (1965). Table of integrals, series and products. Academic Press.
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
require(geoR)

data("Missouri")
data=Missouri
data$V3=log((data$V3))



cc=data$V5
y=data$V3
n=127
datare1=data
coords=datare1[,1:2]
data1=data.frame(coords,y)
data1=data1[cc==0,]
geodata=as.geodata(data1,y.col=3,coords.col=1:2)
v=variog(geodata)
v1=variofit(v)
cov.ini=c(0,2)

est=SAEMSCL(cc,y,cens.type="left",trend="cte",coords=coords,M=15,perc=0.25,MaxIter=5,pc=0.2,
cov.model="exponential",fix.nugget=T,nugget=2,inits.sigmae=cov.ini[2],inits.phi=cov.ini[1],
search=T,lower=0.00001,upper=50)


d1=derivQfun(est)
d1$QI
}

}


\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}
