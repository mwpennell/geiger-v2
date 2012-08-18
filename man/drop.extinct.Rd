\name{drop.extinct}
\alias{drop.extinct}
\alias{drop.random}
\title{ prune specified taxa from a phylogenetic tree }
\description{
  pruning a set of taxa from a tree
}
\usage{
drop.extinct(phy, tol = .Machine$double.eps^0.5) 
drop.random(phy, n)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{tol}{ rounding-error tolerance for taxa that do not reach the present day exactly}
  \item{n}{number of random taxa to prune from the tree}
}

\details{
The functions prune taxa from a tree either at random or based either on a temporal criterion (whether the leaves reach the present within a given \code{tol}). 
The result is a tree that has been pruned based on the given criterion.}

	
\author{ Luke J. Harmon }
\seealso{
	\code{\link[ape]{drop.tip}}
}
\examples{

# Birth-death tree with extinct taxa
p2<-birthdeath.tree(b=0.2, d=0.1, time.stop=30)
plot(p2)

# Previous tree with extinct taxa removed
drop.extinct(p2)->p3
if(!is.null(p3)) plot(p3)

p4<-birthdeath.tree(b=0.2, d=0, taxa.stop=100)
p5<-drop.random(p4, 40)

plot(p5)


}

\keyword{ arith }
