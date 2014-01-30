\name{tips}
\alias{tips}

\title{descendents of a given node in a phylogenetic tree}
\description{
  finding descendants of a node in a tree
}
\usage{
tips(phy, node)
}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{node}{numeric node-identifier (an integer)}
  
}

\value{
The function returns the set of tips subtended by the given \code{node}.
}

\author{ LJ Harmon }

\examples{
geo <- get(data(geospiza))

tips(geo$phy, 18)
}

\keyword{graphs}



