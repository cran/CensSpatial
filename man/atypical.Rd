\name{atypical}
\alias{atypical}

\title{Detection for local influence}
\description{
This function obtain the atypical points detected for the function localinfmeas
}
\usage{
atypical(w)
}

\arguments{
  \item{w}{
An object of the class localinfmes (see the function localinfmeas)
}
}

\value{
\item{RP}{
    Atypical points detected for the response variable perturbation scheme
}


\item{SP}{
    Atypical points detected for the scale matrix perturbation scheme
}

\item{EP}{
    Atypical points detected for the explanatory variable perturbation scheme
}

}
\references{
Cook, R. D. (1986). Assessment of local influence. Journal of the Royal Statistical Society, Series B,, 48, 133-169.

Zhu, H., Lee, S., Wei, B. & Zhou, J. (2001). Case-deletion measures for models with incomplete data. Biometrika,
88, 727-737.
}
\author{

Alejandro Ordonez <<ordonezjosealejandro@gmail.com>>, Victor H. Lachos <<hlachos@ime.unicamp.br>>
and Christian E. Galarza <<cgalarza88@gmail.com>>

Maintainer: Alejandro Ordonez <<ordonezjosealejandro@gmail.com>>
}

\seealso{
\code{\link{localinfmeas}}
}

\examples{
\dontrun{
n<-200 ### sample size for estimation
n1=100 ### number of observation used in the prediction

###simulated coordinates
r1=sample(seq(1,30,length=400),n+n1)
r2=sample(seq(1,30,length=400),n+n1)
coords=cbind(r1,r2)

coords1=coords[1:n,]

cov.ini=c(0.2,0.1)
type="exponential"
xtot=as.matrix(rep(1,(n+n1)))
xobs=xtot[1:n,]
beta=5

###simulated data
obj=rspacens(cov.pars=c(3,.3,0),beta=beta,x=xtot,coords=coords,cens=0.25,n=(n+n1),
n1=n1,cov.model=type,cens.type="left")

data2=obj$datare
cc=obj$cc
y=obj$datare[,3]

##### generating atypical observations###
y[91]=y[91]+4
y[126]=y[126]+4
y[162]=y[162]+4
coords=obj$datare[,1:2]

###initial values###
cov.ini=c(0.2,0.1)

est=SAEMSCL(cc,y,cens.type="left",trend="cte",coords=coords,
M=15,perc=0.25,MaxIter=10,pc=0.2,cov.model=type,
fix.nugget=T,nugget=0,inits.sigmae=cov.ini[1],
inits.phi=cov.ini[2],search=T,lower=0.00001,upper=50)


w=localinfmeas(est,fix.nugget=T,c=3) ## object of class localinfmeas

atypical(w)
}
}

\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}

