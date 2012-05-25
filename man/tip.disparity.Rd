\name{tip.disparity}
\alias{tip.disparity}
\title{ Tip disparity calculation }
\description{
Calculates disparity for every clade in a tree
}
\usage{
tip.disparity(phy, data, data.names=NULL, disp = "avg.sq")
}
\arguments{
  \item{phy}{ object of type 'phylo' }
  \item{data}{Data matrix - either actual or simulated}
  \item{data.names}{Tip names for data vector that match tree species; ignored if data includes names}
  \item{disp}{ Disparity metric }
}
\details{
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
