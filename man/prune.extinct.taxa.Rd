\name{prune.extinct.taxa}
\alias{prune.extinct.taxa}
\alias{prune.random.taxa}
\title{ Prune specified taxa from a phylogenetic tree }
\description{
  Prunes a set of taxa from a tree, either a random fraction of all taxa, or just the taxa whos tip branch 
  terminates before the present day (extinct taxa)
}
\usage{
prune.extinct.taxa(phy, tol = .Machine$double.eps^0.5)
prune.random.taxa(phy, n)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{ Phylogenetic tree }
  \item{tol}{ Tolerance for taxa that do not reach the present day exactly.  This prevents taxa from
  being pruned just due to rounding error}
  \item{n}{Number of random taxa to prune from the tree}
}
\details{
}
\value{
	New tree without these taxa
}
\references{ }
\author{ Luke J. Harmon }
\examples{

# Birth-death tree with extinct taxa
p2<-birthdeath.tree(b=0.2, d=0.1, time.stop=30)
plot(p2)

# Previous tree with extinct taxa removed
prune.extinct.taxa(p2)->p3
if(!is.null(p3)) plot(p3)

p4<-birthdeath.tree(b=0.2, d=0, taxa.stop=100)
p5<-prune.random.taxa(p4, 40)

plot(p5)


}

\keyword{ arith }
