\name{rate.estimate}
\alias{rate.estimate}
\alias{stem.p}
\alias{crown.p}
\alias{crown.limits}
\alias{stem.limits}
\title{ Calculate net diversification rate with confindence limits, and test diversities }
\description{
	rate.estimate Uses Magellon and Sanderson method to calculate net diversification rate for a clade
	given extant diversity and age.
	Associated functions crown.p and stem.p also calculate the probability of obtaining a clade with at least k species given 
	a net diversification rate (r), extinction fraction (e), and time interval.
	Associated functions stem.limits and crown.limits generate confidence limits on extant diversity 
	given a net diversification rate (r), extinction fraction (e), and time interval.
}
\usage{
rate.estimate(time=0, n=0, phy=NULL, epsilon = 0, missing = 0, crown=TRUE, kendall.moran=FALSE)
crown.p(time, r, epsilon, n)
stem.p(time, r, epsilon, n)
crown.limits(r, epsilon, time, prob=c(0.025, 0.975))
stem.limits(r, epsilon, time, prob=c(0.025, 0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{ Time interval; can be a vector }
  \item{n}{Number of extant species}
  \item{phy}{ Phylogenetic tree; can be supplied instead of time and n.  
  			If you're using a tree then crown is automatically true.}
  \item{epsilon}{ Extinction rate as a fraction of speciation rate }
  \item{missing}{ Number of taxa missing from tree }
  \item{crown}{ If true, time is treated as crown age; otherwise, stem age}
  \item{kendall.moran}{ If true, calculates Kendall-Moran estimate of speciation rate; requires a complete phylogenetic tree}
  \item{r}{ Net diversification rate, birth - death }
  \item{prob}{ Range of probabilities for calculating confidence region }
}
\details{
}
\value{
	rate.estimate: Returns net diversification rate r = lambda - mu, and confidence interval
	crown.p and stem.p: Returns the probability of obtaining a clade as big as, 
				or bigger than, size n, given time, r, and epsilon
	stem.limits and crown.limits: Return confidence intervals for clade size given time,
				r, and epsilon
}
\references{ Magallon, S. and M. J. Sanderson. 2000. Absolute diversification rates in angiosperm clades. Evolution 55:1762-1780.}
\author{ Luke J. Harmon and Chad Brock }
\examples{


data(geospiza)
attach(geospiza)

# Assuming no extinction
rate.estimate(phy=geospiza.tree, missing=1)

# Assuming high extinction
rate.estimate(phy=geospiza.tree, epsilon=0.9, missing=1)


}

\keyword{ arith }
