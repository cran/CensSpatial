\name{distmatrix}
\alias{distmatrix}
\title{Distance matrix}
\description{
It computes the euclidean distance matrix for a set of coordinates.
}
\usage{
distmatrix(coords)
}
\arguments{
  \item{coords}{
2D spatial coordinates.
}
}

\value{
\item{dist}{symetric matrix of distances between points.}
}
\references{
Diggle, P. & Ribeiro, P. (2007). Model-Based Geostatistics.
Springer Series in Statistics.
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
n<-200
n1=100

####Simulating spatial coordinates##
r1=sample(seq(1,30,length=400),n+n1)
r2=sample(seq(1,30,length=400),n+n1)
coords=cbind(r1,r2)

H=distmatrix(coords)
}


\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}
