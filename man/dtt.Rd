\name{dtt}
\alias{dtt}
\alias{disparity}
\title{disparity-through-time}
\description{
	calculating and plotting disparity-through-time for a phylogenetic tree and 
	phenotypic data
}
\usage{
disparity(phy=NULL, data, disp = c("avg.sq", "avg.manhattan", "num.states"))
dtt(phy, data, disp=c("avg.sq", "avg.manhattan", "num.states"), nsims=1000, 
	mdi.range=c(0,1), plot=TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{data}{a named vector or matrix of continuous trait values, associated with species in \code{phy}}
  \item{disp}{disparity measure (see \bold{Details})}
  \item{nsims}{number of simulations used to calculate null disparity-through-time plot (see \bold{Details})}
  \item{mdi.range}{time range over which to calculate MDI statistic}
  \item{plot}{whether to plot disparities through time}
}
\details{
		The most complete implementation of \code{dtt} (where \code{nsims} is greater than 0) carries out the entire disparity-through-time (DTT) procedure 
		described in Harmon et al. 2003. The function \code{disparity} simply computes the morphological disparity for a set of species.
		For either function, the disparity measure can be one of the following:
		\itemize{
			\item{\bold{avg.sq }}{is average squared euclidean distance among all pairs of points; this is a good 
					choice if all of the axes are in the same units, or PC axes}
			\item{\bold{avg.manhattan }}{is average Manhattan distance among all pairs of points; this is a good 
					choice if the axes are all in different units, like colors and lengths}
			\item{\bold{num.states }}{is number of unique character states; this is currently the only option for discrete 
					character data}
		}
		
		Note that for \code{mdi.range}, time is relative to a total tree length of 1. The default \code{mdi.range} is the whole tree from 0 (root) to 1 (tips). 
}
\value{
	The function \code{disparity} returns the disparity of the supplied \code{data}. If given a tree, \code{disparity} will return 
	disparities computed for each subtree.  The vector of disparities is indexed based on the numeric node-identifier of the subtending subtree (e.g., the 
	root is indexed by N+1, where N is the number of species in the tree). 
	
	The function \code{dtt} returns elements from the list below:
	\itemize{
		\item{\bold{dtt }}{is average disparity for clades whose stem crosses each time interval in the tree}
		\item{\bold{dtt_sims }}{are disparities at time intervals for each simulated dataset}
		\item{\bold{times }}{are times for each value in disparity-through-time calculation; these are just the branching times of the phylogeny}
		\item{\bold{MDI }}{is the value of the MDI statistic, which is the area between the DTT for the data and the mean of the simulations}
	}	
	If results are plotted, the mean DTT from the simulated datasets appears in dashed line and the empirical DTT in solid line.  
}
\references{
Foote, M. 1997. The evolution of morphological diversity. \emph{ARES} 28:129-152.
		
Harmon, L. J., J. A. Schulte, J. B. Losos, and A. Larson. 2003. Tempo and mode of evolutionary radiation in iguanian lizards. \emph{Science} 301: 961-964.
}

\author{ Luke J. Harmon }
\examples{

geo=get(data(geospiza))

## disparity (not tree-based)
disparity(data=geo$dat) # not tree-based\

## cladewise disparities
disparity(phy=geo$phy, data=geo$dat)

## disparity through time
dttgeo<-dtt(phy=geo$phy, data=geo$dat, nsims=200, plot=TRUE)

}

\keyword{ arith }
