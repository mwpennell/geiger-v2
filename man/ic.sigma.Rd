\name{ic.sigma}
\alias{ic.sigma}
\title{Estimate Evolutionary VCV Matrix}
\description{
	Uses independent contrasts to estimate the evolutionary variance-covariance matrix for 
	n quantitative characters.
}
\usage{
ic.sigma(phy, data, data.names=NULL)
}
\arguments{
  \item{phy}{Phylogenetic tree in 'phylo' format}
  \item{data}{Data matrix}
  \item{data.names}{Tip names for data vector that match tree species; ignored if data includes names}
}
\details{
Data for this function should be a matrix of continuously-valued variables.
}
\value{
	Returns the estimated evolutionary variance-covariance matrix of the variables under a multivariate 
	Brownian motion model.  If you have n characters in your analysis, this will be an nxn matrix.  
	Diagonal elements represent rate estimates for individual characters, while off-diagonal elements
	represent the estimated covariance between two characters.
}
\references{ Revell, L. J., L. J. Harmon, R. B. Langerhans, and J. J. Kolbe. 2007. A phylogenetic 
             approach to determining the importance of constraint on phenotypic evolution in the 
             neotropical lizard, Anolis cristatellus. Evolutionary Ecology Research 9: 261-282. }
\author{ Luke J. Harmon }
\examples{

data(geospiza)
attach(geospiza)

ic.sigma(geospiza.tree, geospiza.data)

}

\keyword{ arith }
