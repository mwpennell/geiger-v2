\name{drop.extinct}
\alias{drop.extinct}
\alias{is.extinct}
\alias{drop.random}
\title{prune specified taxa from a phylogenetic tree}
\description{
  pruning a set of taxa from a tree
}
\usage{
drop.extinct(phy, tol = NULL)
drop.random(phy, n)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{tol}{rounding-error tolerance for taxa that do not reach the present day exactly}
  \item{n}{number of random taxa to prune from the tree}
}

\details{
The functions prune taxa from a tree either at random or based either on a temporal criterion (whether the leaves reach the present within a given \code{tol}). By default, \code{tol = min(phy$edge.length)/100}. The result is a tree that has been pruned based on the given criterion. The function \code{is.extinct} will return a vector of the tip labels of inferred extinct taxa (or NULL is no extinct taxa exist).}


\author{
LJ Harmon, and JW Brown
}
\seealso{
	\code{\link[ape]{drop.tip}}
}
\examples{

# Birth-death tree with extinct taxa
p1 <- sim.bdtree(b=0.2, d=0.1, stop="time", seed=1, t=30)
plot(p1, cex=0.25)

# List extinct taxa
p1.extinct <- is.extinct(p1)

# Previous tree with extinct taxa removed
p2 <- drop.extinct(p1)
plot(p2, cex=0.5)


}

\keyword{graphs}
\keyword{manip}
