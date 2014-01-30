\name{sim.bdtree}
\alias{sim.bdtree}
\title{birth-death tree simulator}
\description{
	simulating phylogenetic trees under a uniform birth-death process
}
\usage{
sim.bdtree(b=1, d=0, stop=c("taxa", "time"), n=100, t=4, seed=0, extinct=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{per-lineage birth (speciation) rate}
  \item{d}{per-lineage death (extinction) rate}
  \item{stop}{stopping criterion}
  \item{n}{maximum number of taxa in simulation}
  \item{t}{maximum time steps of simulation}
  \item{seed}{random number seed (default is to seed based on the clock)}
  \item{extinct}{ whether to return trees where all lineages have gone extinct (see \bold{Details})}  }
\details{
	Starting from a root node - i.e., two living lineages - this function simulates the growth of a
	phylogenetic tree under a uniform, time-homogeneous birth-death process.  This means that every
	lineage has a constant probability of speciating, and a constant probability of going extinct, per
	unit time.  If birth is greater than death, then the number of lineages is expected to grow exponentially.
	If \code{extinct=FALSE}, the function will build trees until one is simulated with at least one surviving lineage. 
}
\value{
	A phylogenetic tree in 'phylo' format is returned.  If death rate is non-zero, then the returned tree will likely 
	include some extinct lineages (terminating before the present day).  See \code{\link[geiger]{drop.extinct}} for 
	a function to remove these lineages.
}
\author{ LJ Harmon and J Weir}
\note{ 
	One note of caution: it is easy to set parameter values that result in tremendously \emph{huge} trees.  If the
	function seems to hang up, this could be the problem.
    
    Other tree simulators are available from the packages \pkg{ape} (\code{\link[ape]{rbdtree}}), \pkg{TreeSim} (\code{\link[TreeSim]{sim.bd.taxa}}, \code{\link[TreeSim]{sim.bd.age}}, and \code{\link[TreeSim]{sim.bd.taxa.age}}), and \pkg{phytools} (\code{\link[phytools]{pbtree}}).
}
\seealso{
\code{\link[geiger]{sim.bd}} for non-phylogenetic simulations; \code{\link[geiger]{drop.extinct}}
}
\examples{

# Pure-birth tree
p1 <- sim.bdtree(b=0.1, d=0, stop="time", t=20)
plot(p1)

# Birth-death tree with extinct taxa
# The extinct flag prevents trees with no survivors

p2 <- sim.bdtree(b=0.2, d=0.05, stop="time", t=20, extinct=FALSE)
plot(p2)

# Previous tree with extinct taxa removed
p3 <- drop.extinct(p2)

}
\keyword{arith}
\keyword{graphs}
