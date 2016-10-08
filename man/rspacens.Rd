\name{rspacens}
\alias{rspacens}

\title{Censored Spatial data simulation

}
\description{

It simulates spatial data with linear structure for one type of censure (left or right).
}
\usage{
rspacens(cov.pars,beta,x=as.matrix(rep(1,n)),coords,kappa=0,cens,n,n1,
cov.model="exponential",cens.type)
}

\arguments{
  \item{cov.pars}{

covariance structure parameters for the errors distribution (\eqn{\phi,  \sigma^2, \tau^2}).
}
  \item{beta}{

linear regression parameters.
}
  \item{x}{

design matrix.
}
  \item{coords}{

coordinates of simulated data.
}
  \item{kappa}{

\eqn{\kappa} parameter used in some covariance structures.
}
  \item{cens}{

percentage of censure in the data (number between 0 and 1).
}
  \item{n}{

number of simulated data used in estimation.
}
  \item{n1}{

number of simulated data used for cross validation (Prediction).
}
  \item{cov.model}{

covariance structure for the data (see \code{cov.spatial} from \code{geoR}).
}
  \item{cens.type}{

 type of censure (\code{"left"} or \code{"right"}).
}
}
\details{
This function analyses prediction in spatial data. It returns a
spatial dataset for estimation (n length) and a spatial dataset (n1 length) used to evaluate the prediction power of a model through
cross validation. The covariance functions used here were provided by cov.spatial from the geoR package.
}

\value{
\item{y}{
    complete simulated data (\eqn{(n+n1)} length).
}

\item{datare}{
    data frame that will be used for the model estimation (coordinates and response).
}

\item{valre}{
    data that will be used for cross validation studies (just response).
}

\item{cc}{
    indicator of censure (1:censored 0:observed).
}

\item{cutof}{
    limit of detection simulated for censure (left: <=cutoff, right: > cutoff).
}

\item{coords1}{
    coordinates of value data.
}



}
\references{

Diggle, P. & Ribeiro, P. (2007). Model-Based Geostatistics. Springer Series in Statistics.

Schelin, L. & Sjostedt-de Luna, S. (2014). Spatial prediction in the presence of left-censoring. Computational Statistics and Data Analysis, 74.


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

n<-200 ### sample size for estimation.
n1=100 ### number of observation used in the prediction.

###simulated coordinates
r1=sample(seq(1,30,length=400),n+n1)
r2=sample(seq(1,30,length=400),n+n1)
coords=cbind(r1,r2)### total coordinates (used in estimation and prediction).
coords1=coords[1:n,]####coordinates used for estimation.

type="matern"### covariance structure.

xtot<-cbind(1,runif((n+n1)),runif((n+n1),2,3))## X matrix for estimation and prediction.
xobs=xtot[1:n,]## X matrix for estimation

obj=rspacens(cov.pars=c(3,.3,0),beta=c(5,3,1),x=xtot,coords=coords,
kappa=1.2,cens=0.25,n=(n+n1),n1=n1,cov.model=type,cens.type="left")

}

\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}

