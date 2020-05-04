\name{predgraphics}
\alias{predgraphics}

\title{Prediction graphics for SAEM Algortihm for censored spatial data.}

\description{
This function provides prediction raster graphics representation and its standard deviation.
}

\usage{
predgraphics(xpred = NULL, grid1, est, points = TRUE,obspoints = 1:sum(est$cc == 0),
colors = terrain.colors(100),sdgraph = TRUE,xlab="X Coord",ylab="Y Coord",
main1="Predicted response", main2="Standard deviation predicted",
xlim=c(min(est$coords[,1]),max(est$coords[,1])),ylim=c(min(est$coords[,2]),
max(est$coords[,2])))
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

\item{xlab}{label for x coordinate of the two plots.}

\item{ylab}{label for y coordinate.}

\item{main1}{an overall title for the prediction plot.}

\item{main2}{an overall title for the standard deviation prediction plot.}

\item{xlim}{x axis limits for the two plots.}

\item{ylim}{y axis limits for the two plots.}

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
\dontshow{
data(depth)
cc=depth$cc[1:70]
y=depth$depth[1:70]
coords=depth[1:70,1:2]

cov.ini=c(1500,30)
est=SAEMSCL(cc,y,cens.type="left",trend="cte",coords=coords,M=10,perc=0.25,
MaxIter=1,pc=0.2,cov.model="gaussian",fix.nugget=FALSE,nugget=10,
inits.sigmae=cov.ini[2],inits.phi=cov.ini[1], search=TRUE,lower=c(0.00001,0.00001),
upper=c(10000,100))


coorgra1=seq(min(coords[,1]),max(coords[,1]),length=30)
coorgra2=seq(min(coords[,2]),max(coords[,2]),length=30)

grid1=expand.grid(x=coorgra1,y=coorgra2)
xpred=rep(1,900)

predgraphics(xpred=xpred,est=est,grid1=grid1,points=TRUE,sdgraph=TRUE)
}

\donttest{
data(depth)
cc=depth$cc
y=depth$depth
coords=depth[,1:2]

cov.ini=c(1500,30)
est=SAEMSCL(cc,y,cens.type="left",trend="cte",coords=coords,M=15,perc=0.25,
MaxIter=100,pc=0.2,cov.model="gaussian",fix.nugget=FALSE,nugget=10,
inits.sigmae=cov.ini[2],inits.phi=cov.ini[1], search=TRUE,lower=c(0.00001,0.00001),
upper=c(10000,100))


coorgra1=seq(min(coords[,1]),max(coords[,1]),length=50)
coorgra2=seq(min(coords[,2]),max(coords[,2]),length=50)

grid1=expand.grid(x=coorgra1,y=coorgra2)
xpred=rep(1,2500)

predgraphics(xpred=xpred,est=est,grid1=grid1,points=TRUE,sdgraph=TRUE)

}
}

\keyword{Spatial}
\keyword{Censored}
\keyword{SAEM}
