\name{area.between.curves}
\alias{area.between.curves}
\title{Area between two curves}
\description{
	Finds the area between two curves (f1 and f2) for a given range on the x-axis.
	One use for this function is to calculate the MDI statistic from Harmon et al. 2003
}
\usage{
area.between.curves(x, f1, f2, xrange = c(0,1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of values for the x axis}
  \item{f1}{Y values for curve 1}
  \item{f2}{Y values for curve 2}
  \item{xrange}{Range of x values over which area is calculated}
}
\details{
	Area is calculated between the two curves f1 and f2 for x values in the given range.
	If f2 is greater than f1, areas are positive; otherwise, areas are negative
}
\value{
	Area between f1 and f2 between 0 and cutoff
}
\references{ Harmon, L. J., J. A. Schulte, J. B. Losos, and A. Larson.  2003.  Tempo and mode of evolutionary radiation in iguanian lizards.  Science 301: 961-964. }
\author{Luke J. Harmon}
\examples{

data(geospiza)
attach(geospiza)

gg<-dtt.full(geospiza.tree, geospiza.data)

# Full area, as in dtt.full function
area.between.curves(gg$times, apply(gg$dtt.sims,1,mean), gg$dtt.data)

# Area for first 2/3, as in Harmon et al. 2003
area.between.curves(gg$times, apply(gg$dtt.sims,1,median), gg$dtt.data, xrange=c(0, 2/3))

# You could also use the mean
area.between.curves(gg$times, apply(gg$dtt.sims,1,mean), gg$dtt.data, xrange=c(0, 2/3))

}

\keyword{arith}
