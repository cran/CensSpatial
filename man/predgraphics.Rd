\name{predgraphics}
\alias{predgraphics}

\title{Prediction graphics for SAEM Algortihm for censored spatial data.

}
\description{
This function provides prediction raster graphics representation and its standard deviation.
}

\usage{
predgraphics(xpred = NULL, grid1, est, points = T,obspoints = 1:length(est$cc == 0),
colors = terrain.colors(100),sdgraph = T,
legend.args.pred = list(text = "Predicted response",side = 4, font = 2, line = 2.5,
cex = 0.8),legend.args.sdpred = list(text = "Standard deviation predicted",side = 4,
font = 2, line = 2.5, cex = 0.8))
}

\arguments{
  \item{xpred}{

x design matrix for the prediction coordinates (must be specified when est$trend="other").
}
  \item{grid1}{

grid with the coordinates of the prediction graphics.
}
  \item{est}{
object of class "SAEMSpatialCens".
}
  \item{points}{
(logical), it indicates if some of the observed points may be plotted in the prediction raster graphic (default, points=\code{TRUE}).
}
  \item{obspoints}{
(vector) if points=\code{TRUE}, it indicates which of the observed (not censored) values may be plotted in the prediction raster graphics.
}
  \item{colors}{

colors pallete used for the graphics (By default \code{terrain.colors(100)}).
}
  \item{sdgraph}{

(logical) it indicates if the standard deviation of the prediction points graphic must be plotted (default sdgraph=\code{TRUE}).
}
  \item{legend.args.pred}{

arguments for the legend of the prediction raster (see \code{plot.raster}).
}
  \item{legend.args.sdpred}{
arguments for the legend of the standard deviation prediction raster (see \code{plot.raster}).
}
}

\value{
in addition to the raster graphics for prediction, the next values are retorned:

\item{datapred}{data.frame with the coordinates and the predicted points used in the prediction raster graphic.}
\item{datasdpred}{data.frame with the coordinates and the standard deviation predicted points used in the standard deviation prediction raster graphic.}


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
\code{\link{SAEMSCL}}
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

type="matern"
xtot=as.matrix(rep(1,(n+n1)))
xobs=xtot[1:n,]
covpars=c(.3,3,0)
beta=5


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

#####predict with a trend in function of its coordinates
xtotcor=as.matrix(cbind(1,coords[,1]))
xcoord=xtotcor[1:n,]


est=SAEMSCL(cc,y,cens.type="left",trend="other",x=xcoord,coords=coords,kappa=1.2,M=15,
perc=0.25,MaxIter=10,pc=0.2,cov.model="exponential", fix.nugget=T,nugget=0,
inits.sigmae=cov.ini[2],inits.phi=cov.ini[1],search=T,lower=0.00001,upper=50)

coorgra1=seq(min(coords[,1]),max(coords[,2]),length=50)
coorgra2=seq(min(coords[,1]),max(coords[,2]),length=50)

grid1=expand.grid(x=coorgra1,y=coorgra2)
xpred=cbind(1,grid1$x)

e=predgraphics(xpred=xpred,grid1=grid1,est,points=T,obspoints=1:50,sdgraph=T)
}
}

\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}
