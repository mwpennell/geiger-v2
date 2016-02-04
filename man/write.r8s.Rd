\name{write.r8s}
\alias{write.r8s}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{
call r8s from geiger}
\description{
call r8s, including a calibration file}
\usage{
write.r8s(phy=NULL, calibrations, base="", blformat=c(lengths="persite", nsites=1, ultrametric="no", round="yes"), divtime=c(method="NPRS", algorithm="POWELL"), describe=c(plot="chrono_description"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogram to turn into a chronogram}
  \item{calibrations}{a set of calibrations}
  \item{base}{idk}
  \item{blformat}{how r8s should interpret your starting tree's branch lengths}
  \item{divtime}{how r8s should make your tree ultrametric}
  \item{describe}{the output you want from r8s}
}
\details{
	This function uses r8s and a calibration to make your tree ultrametric}

\references{
SANDERSON r8s ____NEED TO ADD_____
}

\author{JM Eastman}

\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphs}
\keyword{manip}

