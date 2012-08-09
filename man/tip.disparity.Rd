\name{tip.disparity}
\alias{tip.disparity}
\title{ Tip disparity calculation }
\description{
Calculates disparity for every clade in a tree
}
\usage{
tip.disparity(phy, data, disp = c("avg.sq", "avg.manhattan", "num.states"))
}
\arguments{
  \item{phy}{ object of type 'phylo' }
  \item{data}{Data matrix - either actual or simulated}
  \item{disp}{ Disparity metric }
}
\value{
 	Vector of disparities
}
\references{ Harmon et al. 2003 }
\author{ Luke J. Harmon }
\examples{

data(geospiza)
attach(geospiza)

tip.disparity(geospiza.tree, geospiza.data)

}

\keyword{ arith }
