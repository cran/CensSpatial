\name{depth}
\alias{depth}
\docType{data}
\title{
Depths of a geological horizon.
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Dataset previously analyzed by Dubrule and
Kostov (1986) and De Oliveira (2005).
}
\usage{data("depth")}
\format{
  A data frame with 100 observations on the following 6 variables.
  \describe{
    \item{\code{coord x}}{x coordinate for depth data.}
    \item{\code{coord y}}{y coordinate for depth data.}
    \item{\code{cc}}{indicator of censure (left and right censure).}
    \item{\code{LI}}{lower limit of censure for depth data.}
    \item{\code{LS}}{upper limit of censure for depth data.}
    \item{\code{depth}}{observated depth.}
  }
}

\details{
The observations are placed over a region of about 9 by 5 km and
represent depths of a geological horizon measured at 100 locations where 69 points are fully observed
and 31 points are censored points, these are divided into left- and right- censored points.
The depth data were transformed and their original units remains unknown for confidentiality reasons. For additional details about this dataset we refer to De Oliveira
(2005).

}

\references{
%%  ~~ possibly secondary sources and usages ~~
Dubrule, O. and C. Kostov (1986). An interpolation method taking into account inequality constraints:
I. methodology. Mathematical Geology 18(1), 33-51.

De Oliveira, V. (2005). Bayesian inference and prediction of Gaussian random fields based on censored
data. Journal of Computational and Graphical Statistics 14(1), 95-115.

}
\examples{
data(depth)
summary(depth$depth)
}
\keyword{datasets}
