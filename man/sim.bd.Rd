\name{sim.bd}
\alias{sim.bd}
\title{birth-death population simulator}
\description{
	simulating species richness (or population growth) under a uniform, time-homogeneous 
	birth-death process	
}
\usage{
sim.bd(b=1, d=0, n0=1, times=0:4) 
}
\arguments{
  \item{b}{per-lineage birth (speciation) rate}
  \item{d}{per-lineage death (extinction) rate}
  \item{n0}{number of taxa at starting time zero}
  \item{times}{vector of times where extant species are counted}
}
\details{
	This function simulates species diversification under a uniform birth-death process. This differs 
	from birthdeath.tree in that only the number of species, and not their phylogenetic affinities, are
	stored. This function relates to GEIGER's rate.estimate (and associated functions), which are also non-phylogenetic.
}
\value{
	Population size at each time in TIMES
	}
\references{Yule, G. U. 1924. A mathematical theory of evolution based on the conclusions of 
			Dr. J. C. Willis, FRS. \emph{Philos. Trans. R. Soc. London Ser. B} 213:21-87 }
\author{Richard E. Glor and Luke J. Harmon}
\seealso{
\code{\link[geiger]{bd.ms}}; \code{\link[geiger]{bd.km}}
}

\examples{

pop1=sim.bd(b=0.1, d=0, n0=10, times=1:10)
pop2=sim.bd(b=0, d=0.1, n0=10, times=1:10)
pop3=sim.bd(b=0.1, d=0.1, n0=10, times=1:10)

plot(pop1, type="l", ylim=c(0,max(c(pop1[,"n"], pop2[,"n"], pop3[,"n"]))))
lines(pop2, col="red")
lines(pop3, col="blue")

}
\keyword{arith}
