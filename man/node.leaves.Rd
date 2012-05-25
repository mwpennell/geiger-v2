\name{node.leaves}
\alias{node.leaves}

\title{Returns all descendents of a given node in your tree}
\description{
  Given a node number (as stored in the phylo object's edge matrix), this function returns a list of 
  all descendents of that node.
}
\usage{
node.leaves(phy, node)
}
\arguments{
  \item{phy}{ An object of class phylo}
  \item{node}{Node number}
  
}

\value{
A list of tip names
}

\author{Luke J. Harmon}

\examples{
data(geospiza)
attach(geospiza)

node.leaves(geospiza.tree, 18)
}

\keyword{manip}



