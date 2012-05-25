\name{dtt}
\alias{dtt}
\alias{disp.calc}
\alias{dtt.full}
\title{Disparity-through-time}
\description{
	Functions for calculating and plotting disparity-through-time for a phylogenetic tree and 
	phenotypic data.
}
\usage{
disp.calc(data, disp = "avg.sq")
dtt(phy, data, data.names=NULL, disp = "avg.sq")
dtt.full(phy, data, data.names=NULL, disp="avg.sq", nsims=1000, mdi.range=c(0,1))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data matrix - either actual or simulated}
  \item{data.names}{Tip names for data vector that match tree species; ignored if data includes names}
  \item{phy}{Phylogenetic tree in 'phylo' format}
  \item{disp}{Disparity measure: see below for currently implemented options.}
  \item{nsims}{Number of simulations used to calculate null dtt plot}
  \item{mdi.range}{Time range over which to calculate MDI statistic (area between curves).  Time is 
  					relative to the total tree length of 1; default is the whole tree from 0 (root) to 1 (tips)}
}
\details{
		The most complete function, dtt.plot, carries out the entire disparity-through-time procedure 
		described in Harmon et al. 2003.  Other functions are for various parts of this, as follows:
		disp.calc: Calculate morphological disparity for a set of species.
		Disparity measure can be one of the following:
		avg.sq: Average squared euclidean distance among all pairs of points.  This is a good 
					choice if all of the axes are in the same units, or PC axes.
		avg.manhattan: Average Manhattan distance among all pairs of points.  This is a good 
					choice if the axes are all in different units, like colors and lengths.
		nb.states: Number of unique character states; this is the only option for discrete 
					character data, for now.
		
		dtt: Evaluates disparity-through-time for either a single data set or multiple 
				data sets (for example, from simulations).
		dtt.full: Calculates dtt for the actual data, runs simulations using univariate or multivariate
				Brownian motion, calculates dtt for all simulations and finds the mean, and plots dtt for 
				both data and simulations on the output device.

}
\value{
	disp.calc: Disparity of the supplied data
	dtt: Average disparity for clades whose stem crosses each time interval in the tree.
	dtt.full: A list with the following items:
	dtt.data: DTT for the data
	dtt.sims: DTT for each simulated data set
	times: Times for each value in the dtt plot; this is just the branching times of the phylogeny
	MDI: Value of the MDI statistic, which is the area between the DTT plot for the 
					data and the mean of the simulations
`	Plot: Creates a DTT plot	
	
}
\references{Foote, M. 1997. The evolution of morphological diversity. Annual Review of Ecology and 
			Systematics 28:129-152.
			Harmon, L. J., J. A. Schulte, J. B. Losos, and A. Larson. 2003. Tempo and mode of 
			evolutionary radiation in iguanian lizards. Science 301: 961-964.}
\author{ Luke J. Harmon }
\examples{

data(geospiza)
attach(geospiza)

disp.calc(geospiza.data)
disp.data<-dtt(geospiza.tree, geospiza.data)
full.output<-dtt.full(geospiza.tree, geospiza.data)

}

\keyword{ arith }
