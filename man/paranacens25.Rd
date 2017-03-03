\name{paranacens25}
\alias{paranacens25}
\docType{data}
\title{
Artificially censored Rainfall Data from Parana State, Brasil
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Artificially censored data-set used by Diggle and Ribeiro (2001) to illustrate the methods discussed in the paper. The censure for these data was generated following the methods of Schelin et al. (2014).
}
\usage{data("paranacens25")}
\format{
  A data frame with 143 observations on the following 6 variables.
  \describe{
    \item{\code{east}}{x coordinate for depth data.}
    \item{\code{north}}{y coordinate for depth data.}
    \item{\code{y1}}{indicator of censure (left and right censure).}
    \item{\code{cc}}{lower limit of censure for depth data.}
  }
}


\references{
%%  ~~ possibly secondary sources and usages ~~
Schelin, L. & Sjostedt-de Luna, S. (2014). Spatial prediction in the presence of left-censoring. Computational Statistics and Data Analysis, 74.

Diggle, P.J. \& Ribeiro Jr, P.J. (2002) Bayesian inference in Gaussian model-based geostatistics. Geographical and Environmental Modelling, Vol. 6, No. 2, 129-146.

}
\examples{
data(paranacens25)
summary(paranacens25$y1)
}
\keyword{datasets}
