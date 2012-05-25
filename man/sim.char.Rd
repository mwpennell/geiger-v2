\name{sim.char}
\alias{sim.char}
\title{ Simulate character evolution }
\description{
  Simulates evolution of discrete or continuous characters on a phylogenetic tree
}
\usage{
sim.char(phy, model.matrix, nsims = 1, model = "brownian", root.state = 1)
}
\arguments{
  \item{phy}{ Phylogenetic tree }
  \item{model.matrix}{ Matrix describing model: either vcv matrix or q matrix }
  \item{nsims}{ Number of simulations to run }
  \item{model}{ "brownian", "speciational", or "discrete" }
  \item{root.state}{ Starting state at root, only needed for discrete model. }
}
\details{
	This function simulates either discrete or continuous data on a phylogenetic tree. The model variable 
	determines the type of simulation to be run. There are three options: "discrete," which evolves
	characters under a continuous time Markov model, and two continuous models, "brownian" and "speciational."
	Brownian is a constant rate Brownian motion model, while speciational is a Brownian model on a tree
	where all branches have the same length.  The model.matrix parameter gives the structure of the model, 
	and should be either a transition matrix, Q, for model = "discrete," or a trait variance-covariance 
	matrix for model = "brownian" or "speciational."  For discrete models, one character at a time is 
	simulated, and the size of the model.matrix gives the number of character states; for continuous models, 
	multivariate characters can be simulated, with their evolution goverened by a covariance matrix specified 
	in the model.matrix.  "root.state" is only needed for the discrete model.

}
\value{
 	Array of simulated data, either two or three-dimensional. 
 	The first dimension is the number of taxa, the second the number of characters, and the third the 
 	number of simulated data sets.
}
\references{ }
\author{ Luke J, Harmon }
\examples{
data(geospiza)
attach(geospiza) 


s<-ic.sigma(geospiza.tree, geospiza.data)
sims<-sim.char(geospiza.tree, s, 100)

# Discrete character evolution
q<-list(rbind(c(-.5, .5), c(.5, -.5)))
sims<-sim.char(geospiza.tree, q, model="discrete", n=100)
  
}

\keyword{ arith }
