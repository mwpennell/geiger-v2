\name{sim.char}
\alias{sim.char}
\title{ simulate character evolution }
\description{
  simulating evolution of discrete or continuous characters on a phylogenetic tree
}
\usage{
sim.char(phy, model.matrix, nsims = 1, model = c("brownian", "speciational", "discrete"), root = 1)
}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{model.matrix}{matrix describing model (either vcv matrix or Q matrix) }
  \item{nsims}{ number of simulations to run }
  \item{model}{a model from which to simulate data}
  \item{root}{starting state (value) at root}
}
\details{
	This function simulates either discrete or continuous data on a phylogenetic tree. The model variable 
	determines the type of simulation to be run. There are three options: \code{discrete}, which evolves
	characters under a continuous time Markov model, and two continuous models, \code{brownian} and \code{speciational}.
	\code{brownian} is a constant rate Brownian-motion model, while \code{speciational} is a Brownian model on a tree
	where all branches have the same length.  The \code{model.matrix} parameter gives the structure of the model, 
	and should be either a transition matrix, Q, for the \code{discrete} model, or a trait variance-covariance 
	matrix for \code{brownian} or \code{speciational} models.  For discrete models, multiple characters may be simulated 
	if \code{model.matrix} is given as a list of Q matrices (see \bold{Examples}). For continuous models, multivariate characters can be simulated, 
	with their evolution goverened by a covariance matrix specified in the \code{model.matrix}. 

}
\value{
 	An array of simulated data, either two or three-dimensional, is returned. 
 	The first dimension is the number of taxa, the second the number of characters, and the third the 
 	number of simulated data sets.
}
\author{ Luke J, Harmon }
\examples{
geo=get(data(geospiza))


## Continuous character -- univariate
usims<-sim.char(geo$phy, 0.02, 100)

## Continuous character -- multivariate
s<-ic.sigma(geo$phy, geo$dat)
csims<-sim.char(geo$phy, s, 100)

## Discrete character -- univariate
q<-list(rbind(c(-.5, .5), c(.5, -.5)))
dsims<-sim.char(geo$phy, q, model="discrete", n=10)

## Discrete character -- multivariate
qq<-list(rbind(c(-.5, .5), c(.5, -.5)), rbind(c(-.05, .05), c(.05, -.05)))
msims<-sim.char(geo$phy, qq, model="discrete", n=10)

  
}

\keyword{ arith }
