\name{Missouri}
\alias{Missouri}
\docType{data}
\title{
TCDD concentrations in Missouri (1971).
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Contents the data of TCDD concentrations used for Zirschky et al. in his geostatistical analysis of Hazardous waste data in Missouri.

}
\usage{data("Missouri")}
\format{
  A data frame with 127 observations on the following 5 variables.
  \describe{
    \item{\code{V1}}{x coordinate of start of each transect (ft).}
    \item{\code{V2}}{y coordinate of start of each transect (ft).}
    \item{\code{V3}}{TCDD Concentrations (mg/m^3).}
    \item{\code{V4}}{transect length (ft).}
    \item{\code{V5}}{indicator of censure (left censure in all data).}
  }
}

\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
The data was collected in November 1983 by U.S. EPA in several areas of a highway from Missouri. Only the locations used in the geostatistical analysis  by the autors are showed.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Zirschky, J. H. & Harris, D. J. (1986). Geostatistical analysis of hazardous waste site data. Journal of Environmental
Engineering, 112(4), 770-784.
}
\examples{
data(Missouri)
summary(Missouri$V3)
}
\keyword{datasets}
