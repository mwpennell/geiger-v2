\name{dtt}
\alias{dtt}
\alias{disparity}
\title{disparity-through-time}
\description{
	calculating and plotting disparity-through-time for a phylogenetic tree and 
	phenotypic data
}
\usage{
disparity(phy=NULL, data, index = c("avg.sq", "avg.manhattan", "num.states"))
dtt(phy, data, index=c("avg.sq", "avg.manhattan", "num.states"),
    mdi.range=c(0,1), nsim=0, CI=0.95, plot=TRUE, calculateMDIp=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{data}{a named vector or matrix of continuous trait values, associated with species in \code{phy}}
  \item{index}{disparity index to use (see \bold{Details})}
  \item{mdi.range}{time range over which to calculate MDI statistic}
  \item{nsim}{number of simulations used to calculate null disparity-through-time plot (see \bold{Details})}
  \item{CI}{confidence level from which to plot simulated disparities}
  \item{plot}{whether to plot disparities through time}
  \item{calculateMDIp}{calculate p-value for MDI compared to null; only valid if nsim is not zero}
}
\details{
		The most complete implementation of \code{dtt} (where \code{nsim} is greater than 0) carries out the entire disparity-through-time (DTT) procedure described in Harmon et al. 2003, where simulated data are used to construct a null DTT distribution. The function \code{disparity} simply computes the morphological disparity for a set of species. Note that for \code{mdi.range}, time is relative to a total tree length of 1. The default \code{mdi.range} is the entire temporal span of the tree, from 0 (root) to 1 (tips). 
		
        For either function, the disparity index can be one of the following:
		\itemize{
			\item{\bold{avg.sq }}{is average squared Euclidean distance among all pairs of points; this is the most common distance metric for disparity in macroevolution}
			\item{\bold{avg.manhattan }}{is average Manhattan distance among all pairs of points}
			\item{\bold{num.states }}{is number of unique character states; this is currently the only option for discrete 
					character data}
		}
}
\value{
	The function \code{disparity} returns the disparity of the supplied \code{data}. If given a tree, \code{disparity} will return 
	disparities computed for each subtree.  The vector of disparities is indexed based on the numeric node-identifier of the subtending subtree (e.g., the 
	root is indexed by N+1, where N is the number of species in the tree; see \code{\link[ape]{read.tree}}). 
	
	The function \code{dtt} returns elements from the list below:
	\itemize{
		\item{\bold{dtt }}{is average disparity for clades whose stem crosses each time interval in the tree}
        \item{\bold{times }}{are times for each value in disparity-through-time calculation; these are just the branching times of the phylogeny}
		\item{\bold{sim }}{are disparities at time intervals for each simulated dataset}
		\item{\bold{MDI }}{is the value of the MDI statistic, which is the area between the DTT for the data and the median of the simulations}
	}	
	If results are plotted, the mean DTT from the simulated datasets appears in dashed line and the empirical DTT in solid line.  
}
\references{
Foote M. 1997. The evolution of morphological diversity. \emph{ARES} 28:129-152.
		
Harmon LJ, JA Schulte, JB Losos, and A Larson. 2003. Tempo and mode of evolutionary radiation in iguanian lizards. \emph{Science} 301:961-964.

Slater GJ, SA Price, F Santini, and MA Alfaro. 2010. Diversity vs disparity and the evolution of modern cetaceans. \emph{PRSB} 277:3097 -3104.
}

\author{ LJ Harmon and GJ Slater}
\examples{
\dontrun{
geo=get(data(geospiza))

## disparity -- not tree-based
disparity(data=geo$dat) # not tree-based

## cladewise disparities
disparity(phy=geo$phy, data=geo$dat)

## disparity through time of culmen length
dttcul<-dtt(phy=geo$phy, data=geo$dat[,"culmenL"], nsim=100, plot=TRUE)

## disparity through time of entire dataset -- without simulated data
dttgeo<-dtt(phy=geo$phy, data=geo$dat, nsim=0, plot=TRUE)
}
}

\keyword{ arith }
