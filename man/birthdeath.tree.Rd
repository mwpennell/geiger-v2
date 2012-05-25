\name{birthdeath.tree}
\alias{birthdeath.tree}
\title{Birth-death tree simulator}
\description{
	Simulates phylogenetic trees under a uniform birth-death process.
}
\usage{
birthdeath.tree(b, d, time.stop = 0, taxa.stop = 0, seed = 0, print.seed=FALSE, return.all.extinct=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{Per-lineage birth (speciation) rate}
  \item{d}{Per-lineage death (extinction) rate}
  \item{time.stop}{Stopping time}
  \item{taxa.stop}{Maximum number of taxa}
  \item{seed}{Random number seed; if seed=0 (default) then random number 
              generator is seeded based on the clock}
  \item{print.seed}{If T, prints out the random number seed associated with each tree} 
  \item{ return.all.extinct }{ Return trees where all lineages have gone extinct?  Otherwise, the function will 
  								build trees until it gets one with at least one surviving lineage. }
}
\details{
	Starting from a root node - i.e., two living lineages - this function simulates the growth of a
	phylogenetic tree under a uniform, time-homogeneous birth-death process.  This means that every
	lineage has a constant probability of speciating, and a constant probability of going extinct, per
	unit time.  If birth is greater than death, then the number of lineages is expected to grow exponentially.

}
\value{
	Phylogenetic tree in ape format.  If death rate is non-zero, then the returned tree will likely 
	include some extinct lineages (terminating before the present day).  The GEIGER function 
	prune.extinct.taxa can remove these lineages.
}
\references{Geiger}
\author{Luke J. Harmon and Jason Weir}
\note{ 
	One note of caution: it is easy to set parameter values that result in tremendously HUGE trees.  If the
	function seems to hang up, this could be the problem.
}
\seealso{
BDsim for non-phylogenetic simulations; prune.extinct.taxa
}
\examples{

# Pure-birth tree
p1<-birthdeath.tree(b=0.1, d=0, time.stop=20)
plot(p1)

# Birth-death tree with extinct taxa
# The return.all.extinct flag prevents trees with no survivors

p2<-birthdeath.tree(b=0.2, d=0.05, time.stop=20, return.all.extinct=FALSE)
plot(p2)

# Previous tree with extinct taxa removed
prune.extinct.taxa(p2)->p3

}
\keyword{arith}
\keyword{graphs}
