\name{getAncStates}
\alias{getAncStates}
\title{Estimate ancestral character states for continuous characters}
\description{
	Estimates ancestral character states for continuous characters under 
	a Brownian motion model
}
\usage{
getAncStates(x, phy) 
}
\arguments{
  \item{x}{Data vector}
  \item{phy}{Phylogenetic tree in phylo format}
}
\details{
	This function uses an algorithm to calculate Maximum-likelihood estimates of 
	ancestral character states for continuous characters.
	
	The ape function ace carries out a similar function but finds estimates
	by actually maximizing the likelihood across all ancestors simultaneously.
	
	I have implemented an algorithm that is faster and sometimes more reliable
	especially for large trees.
	}
\value{
	Vector of ancestral character states at each node.
	Nodes are numbered following the order in phy edge
	}
\references{ }
\author{Luke J. Harmon}
\seealso{
ace
}
\examples{

data(geospiza)
attach(geospiza)

tt<-drop.tip(geospiza.tree, "olivacea")
getAncStates(geospiza.data[,1], tt)
 

}
\keyword{arith}
