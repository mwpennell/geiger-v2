\name{bd.ms}
\alias{bd.ms}
\alias{bd.km}
\alias{stem.p}
\alias{crown.p}
\alias{crown.limits}
\alias{stem.limits}
\title{ estimate net diversification rate }
\description{
estimating net diversification rate with confidence limits and testing diversities
}
\usage{
bd.ms(phy=NULL, time, n, missing = 0, crown=TRUE, epsilon = 0)
bd.km(phy=NULL, time, n, missing = 0, crown=TRUE)
crown.p(phy=NULL, time, n, r, epsilon)
stem.p(phy=NULL, time, n, r, epsilon)
crown.limits(time, r, epsilon, CI=0.95)
stem.limits(time, r, epsilon, CI=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo' (see \bold{Details})}
  \item{time}{ time interval (can be a vector) }
  \item{n}{number of extant species}
  \item{r}{net diversification rate, birth - death }
  \item{epsilon}{ extinction rate as a fraction of speciation rate }
  \item{missing}{ number of taxa missing from tree }
  \item{crown}{whether time is treated as crown age (or otherwise as stem age)}
  \item{CI}{confidence level for estimated parameters}
}
\details{	
	\code{bd.ms} uses the Magallon and Sanderson (2000) method to calculate net diversification rate for a clade given extant diversity and age. \code{bd.km} 
	computes the Kendall-Moran estimate of speciation rate, which assumes a complete phylogenetic tree.


	Associated functions \code{crown.p} and \code{stem.p} also calculate the probability of obtaining a clade with at least \code{n} species given 
	a net diversification rate (\code{r}), extinction fraction (\code{epsilon}), and \code{time} interval. Associated functions \code{stem.limits} and 
	\code{crown.limits} generate confidence limits on extant diversity given a net diversification rate (\code{r}), extinction fraction (\code{epsilon}), 
	and \code{time} interval.
	
	Where a function calls for a \code{time} and an \code{n} element, a tree may be given instead (as argument \code{phy}). The argument \code{n} is taken from the number of tips in the tree. The method 
	will attempt to discern whether the model should be fitted assuming \code{crown} or \code{stem}.  If the tree has a non-NULL \code{phy$root.edge} element, the length will be assumed for the stem and \code{crown} is assumed to be \code{FALSE} (see also \code{\link[ape]{read.tree}}). 
}

\value{
	\itemize{
		\item{\bold{bd.ms} }{returns net diversification rate (\code{r} = \code{lambda} - \code{mu})}
		\item{\bold{bd.km} }{returns speciation rate assuming a completely sampled tree}
		\item{\bold{crown.p} and \bold{stem.p} }{return the probability of obtaining a clade as big as (or bigger than) size \code{n}, given 
		\code{time}, \code{r}, and \code{epsilon}}
		\item{\bold{crown.limits} and \bold{stem.limits} }{return lower (\code{lb}) and upper (\code{ub}) confidence intervals for clade size given \code{time}, 
		\code{r}, and \code{epsilon}}	
	}
}

\references{ 
Magallon S and MJ Sanderson. 2000. Absolute diversification rates in angiosperm clades. \emph{Evolution} 55:1762-1780.
}

\author{ LJ Harmon and C Brock }

\examples{


geo=get(data(geospiza))

# Assuming no extinction
bd.ms(phy=geo$phy, missing=1, epsilon=0)

# Assuming high extinction
bd.ms(phy=geo$phy, missing=1, epsilon=0.9)


}

\keyword{ arith }
