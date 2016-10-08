\name{derivcormatrix}
\alias{derivcormatrix}

\title{First and second derivates of some correlation matrix
}
\description{
It computes the matrix of first and second derivates for the
exponential, gaussian, matern, spherical, powered exponential and
Cauchy correlation matrix.
}
\usage{
derivcormatrix(coords, phi, kappa = 0, cov.model = "exponential")
}
\arguments{
  \item{coords}{
2D spatial coordinates.
}
  \item{phi}{
parameter for the matern, powered exponential and cauchy functions.
}
  \item{kappa}{
parameter for all correlation functions.
}
  \item{cov.model}{
parameter correlation funtion to calculates the derivates in this case 6 functions are avalible
"exponential", "gaussian", "matern", "spherical", "powered.exponential", "cauchy".
}
}
\details{
The correlations functions used to calculate the derivates
from this 6 functions are based in the functions
by the package \code{geoR} (see \code{cov.spatial}).
}
\value{
\item{H}{distance matrix.}
\item{devR1}{first derivate of the correlation matrix.}
\item{devR2}{second derivate of the correlation matrix.}
}

\seealso{
\code{\link{SAEMSCL}}
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


\examples{
n<-200
n1=100
r1=sample(seq(1,30,length=400),n+n1)
r2=sample(seq(1,30,length=400),n+n1)
coords=cbind(r1,r2)

s=derivcormatrix(coords=coords,phi=2,kappa=2,cov.model="exponential")

}
\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}
