\name{localinfmeas}
\alias{localinfmeas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Local influence measures.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
It computes some measures and plots to asses the local influence of outliers in the SAEM spatial estimation
for censored spatial observations, for six types of covariance functions (est$type): "exponential", "matern", "gauss", "spherical","powered.exponential" or "stable" and "cauchy".
}
\usage{
localinfmeas(est, fix.nugget = TRUE, diag.plot = TRUE, type.plot = "all", c = 3)
}

\arguments{
  \item{est}{
object of the class "SAEMSpatialCens". See \code{SAEMSCL} function.
}
  \item{fix.nugget}{

(logical) it indicates if the \eqn{\tau^2} parameter must be fixed.
}
  \item{diag.plot}{

(logical) it indicates if diagnostic plots must be showed.
}
  \item{type.plot}{
%%     ~~Describe \code{type.plot} here~~
type of plot (all: all graphics, rp: response perturbation,smp: scale matrix perturbation,
evp: explanatory variable perturbation).
}
  \item{c}{
%%     ~~Describe \code{c} here~~
constant used for fixing the limit of detection (benchmark value).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
this function uses the Maximum likelihood expectation (MLE) under three perturbation schemes
(response, scale matrix and explanatory variables perturbations) to detect the influence of outliers in the SAEM estimation procedure.
}

\value{
in addition to the diagnostic graphics (response, scale matrix and explanatory variable schemes, respectively), the function returns the next values.

\item{Qwrp}{negative \eqn{Qw0} matrix under the response perturbation scheme.}
\item{Qwsmp}{negative \eqn{Qw0} matrix under the scale matrix perturbation scheme.}
\item{Qwevp}{negative \eqn{Qw0} matrix under the explanatory variable perturbation scheme.}
\item{respper}{data.frame containing an indicator of the presence of atypical values
and the \eqn{M(0)} values for the response perturbation scheme.}
\item{smper}{data.frame containing an indicator of the presence of atypical values
and the \eqn{M(0)} values for the scale matrix perturbation scheme.}
\item{expvper}{data.frame containing an indicator of the presence of atypical values
and the \eqn{M(0)} values for the explanatory variable perturbation scheme.}
\item{limrp}{limit of detection for outliers for the response perturbation scheme.}
\item{limsmp}{limit of detection for outliers for the scale matrix perturbation scheme.}
\item{limevp}{limit of detection for outliers for the explanatory variable perturbation scheme.}

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
\code{\link{SAEMSCL}}
}

\examples{
\dontshow{
require(geoR)

data("Missouri")
data=Missouri[1:70,]
data$V3=log((data$V3))
cc=data$V5
y=data$V3
n=127
k=1
datare1=data
coords=datare1[,1:2]
data1=data.frame(coords,y)
data1=data1[cc==0,]
geodata=as.geodata(data1,y.col=3,coords.col=1:2)
v=variog(geodata)
v1=variofit(v)
cov.ini=c(0,2)
est=SAEMSCL(cc,y,cens.type="left",trend="cte",coords=coords,M=15,perc=0.25,
MaxIter=1,pc=0.2,cov.model="exponential",fix.nugget=TRUE,nugget=2,
inits.sigmae=cov.ini[2],inits.phi=cov.ini[1], search=TRUE,lower=0.00001,upper=100)


w=localinfmeas(est,fix.nugget=TRUE,c=3)

res=w$respper
res[res[,1]=="atypical obs",]

sm=w$smper
sm[sm[,1]=="atypical obs",]

ev=w$expvper
ev[ev[,1]=="atypical obs",]
}

\donttest{

require(geoR)

data("Missouri")
data=Missouri
data$V3=log((data$V3))
cc=data$V5
y=data$V3
n=127
k=1
datare1=data
coords=datare1[,1:2]
data1=data.frame(coords,y)
data1=data1[cc==0,]
geodata=as.geodata(data1,y.col=3,coords.col=1:2)
v=variog(geodata)
v1=variofit(v)
cov.ini=c(0,2)
est=SAEMSCL(cc,y,cens.type="left",trend="cte",coords=coords,M=15,perc=0.25,
MaxIter=5,pc=0.2,cov.model="exponential",fix.nugget=TRUE,nugget=2,
inits.sigmae=cov.ini[2],inits.phi=cov.ini[1], search=TRUE,lower=0.00001,upper=100)


w=localinfmeas(est,fix.nugget=TRUE,c=3)

res=w$respper
res[res[,1]=="atypical obs",]

sm=w$smper
sm[sm[,1]=="atypical obs",]

ev=w$expvper
ev[ev[,1]=="atypical obs",]


##############ANOTHER EXAMPLE#########

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

est=SAEMSCL(cc,y,cens.type="left",trend="cte",coords=coords,M=15,perc=0.25,
MaxIter=10,pc=0.2,cov.model=type,fix.nugget=TRUE,nugget=0,inits.sigmae=cov.ini[1],
inits.phi=cov.ini[2],search=TRUE,lower=0.00001,upper=50)


w=localinfmeas(est,fix.nugget=TRUE,c=3)

res=w$respper
res[res[,1]=="atypical obs",]

sm=w$smper
sm[sm[,1]=="atypical obs",]

ev=w$expvper
ev[ev[,1]=="atypical obs",]

}
}


\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}
