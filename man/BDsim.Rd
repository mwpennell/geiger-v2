\name{BDsim}
\alias{BDsim}
\title{Birth-death population simulator}
\description{
	Simulates species richness (or population growth) under a uniform, time-homogeneous 
	birth-death process.	
}
\usage{
BDsim(nStart, b, d, times) 
}
\arguments{
  \item{nStart}{Number of taxa at starting time zero}
  \item{b}{Per-lineage birth (speciation) rate}
  \item{d}{Per-lineage death (extinction) rate}
  \item{times}{Vector of times where extant species are counted}
}
\details{
	This function simulates species diversification under a uniform birth-death process.  This differs 
	from birthdeath.tree in that only the number of species, and not their phylogenetic affinities, are
	stored.  This function relates to GEIGER's rate.estimate (and associated functions), which are also non-phylogenetic.
}
\value{
	Population size at each time in TIMES
	}
\references{Yule, G. U. 1924. A mathematical theory of evolution based on the conclusions of 
			Dr. J. C. Willis, FRS. Philos. Trans. R. Soc. London Ser. B 213: 21-87 }
\author{Richard E. Glor and Luke J. Harmon}
\seealso{
rate.estimate, stem.p, crown.p, stem.limits, crown.limits
}
\examples{

pop1=BDsim(nStart=10, b=0.1, d=0, times=1:10)
pop2=BDsim(nStart=10, b=0, d=0.1, times=1:10)
pop3=BDsim(nStart=10, b=0.1, d=0.1, times=1:10)

plot(pop1, type="l", ylim=c(0,20))
lines(pop2, col="red")
lines(pop3, col="blue")

}
\keyword{arith}
