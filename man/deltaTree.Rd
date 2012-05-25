\name{deltaTree}
\alias{deltaTree}
\alias{lambdaTree}
\alias{kappaTree}
\alias{ouTree}
\alias{tworateTree}
\alias{linearchangeTree}
\alias{rescaleTree}
\alias{exponentialchangeTree}
\alias{speciationalTree}
\title{Tree transformations}
\description{
	Apply various transformation to the branches of a phylogenetic tree.
	}
\usage{
deltaTree(phy, delta, rescale = T)
lambdaTree(phy, lambda)
kappaTree(phy, kappa)
ouTree(phy, alpha)
tworateTree(phy, breakPoint, endRate) 
linearchangeTree(phy, endRate=NULL, slope=NULL)
exponentialchangeTree(phy, endRate=NULL, a=NULL)
speciationalTree(phy)
rescaleTree(phy, totalDepth)
}
\arguments{
  \item{phy}{ an object of class phylo }
  \item{delta}{Delta value}
  \item{rescale}{if TRUE, rescale tree so that the total depth (root-to-tip distance) remains constant; this facilitates comparisons of rates among different delta parameter values}
  \item{lambda}{Lambda value}
  \item{kappa}{ Kappa value }
  \item{alpha}{Value of Ornstein-Uhlenbeck constraint parameter}
  \item{breakPoint}{Applies a new rate of evolution after the time specified as the breakpoint}
  \item{endRate}{Rate of change at the present day, relative to the initial (root) rate}
  \item{slope}{Slope of the relationship between rate and time in the linearchange model}
  \item{a}{Exponent of the relationship between rate and time in the exponentialchange model}
  \item{totalDepth}{Desired value for root to tip distance of the tree.} 
}
\details{
	All of these functions take a tree and various parameters as inputs, and output a new tree structure 
	that has been transformed in the appropriate manner.  These are meant to correspond with changing the model of phenotypic evolution 
	for discrete or continuous characters.  To use the modified model, one then evolves characters on the transformed tree under a 
	Brownian motion model.  Another use for these functions is in finding the particular parameter values that maximize the likelihood
	of the given data; they can then be compared to Brownian motion using a likelihood-ratio test, or to any other model using AIC.
	deltaTree: A transformation suggested by Pagel as a test for a slow-down or speed-up in the rate of character
					evolution through time.  Value of delta > 1 disproporionately increase the length of external nodes 
					(speed-up), where as values of delta < 1 disproportionately increase the length of internal nodes (slow-down).
					Delta = 1 is a Brownian motion model (so the tree is returned unchanged).
	lambdaTree: Another Pagel transformation. This one manipulates the tree as a test of phylogenetic signal.  
				Values of lambda should always less than one.  Each internal edge is multiplied by lambda, while branches
				leading to tips are unchanged; this has the effect of reducing or eliminating phylogenetic signal.	Lambda = 1 
				is a Brownian motion model (so the tree is returned unchanged).
	kappaTree: A third Pagel transformation. This one manipulates the tree as a test of "speciational" models.  
				Each branch length in the tree is raised to the power kappa. Kappa = 1 is a Brownian motion model (so the tree is returned unchanged), 
				while kappa = 0 is a speciational model where all branch lengths are equal.  Kappa > 1 or < 0 does not make sense to me, so I think
				the values should be restricted to between zero and one.
	ouTree: A transformation that corresponds to evolution under an Ornstein-Uhlenbeck model.  This model can be thought of as a 
						random-walk model with a central tendency, so that phynotypes tend to evolve towards one "optimal" value.  Evolving characters
						on the tree that's returned from this function mimics evolving them under the OU model in the case where the optimal value
						is exactly equal to the ancestral state for the character.  If these two values are NOT equal, then simply transforming the tree
						is not enough to run simulations - and I haven't implemented any other way of doing this.
	tworateTree: A transformation that effectively changes the rate of evolution at some point in time, specified by breakPoint.
						The rate after the breakpoint is supplied as endRate, and only needs to be specified as a ratio relative to the initial rate.  
						Thus, if endRate is greater than 1, evolution speeds up (i.e. all branches towards the tips are made longer; 
						if endrate is less than 1, branches towards the tips are made shorter. If endRate = 1, the model is a constant-
						rate Brownian motion model (so the tree is returned unchanged).
	linearchangeTree: A transformation that changes the rate of evolution linearly through time.
						The rate at the present day is supplied either as endRate, a ratio relative to the initial rate
						at the root of the tree, or as the slope of the relationship between rate and time.  Rates change linearly: 
						r(t) = ro + slope * t, where ro is the inital rate and a is the slope determined using endRate. 
						If endRate is greater than 1 or slope > 0, evolution gradually speeds up; if endrate is less than 1 or slope < 0, 
						evolution gradually slows down. If endRate = 1, the model is a constant-
						rate Brownian motion model (so the tree is returned unchanged).
	exponentialchangeTree: A transformation that changes the rate of evolution exponentially through time.
						The rate at the present day is supplied as endRate or a, in exactly the same way as in linearchangeTree.
						The difference for this function is that rates grow or decay exponentially rather than linearly.
						The equation is r(t) = ro * exp(a * t), where ro is the inital rate and a either given or is determined using endRate.
						One of the two (a or endRate) must be specified, or the function returns an error.
	speciationalTree: Transforms all the branch lengths in the tree to a length of one; sometimes this is called a 'speciational' or
						'punctuated' model of evolution.
	rescaleTree: Rescales a tree so that the total depth of the tree (i.e. the distance from the root to the tips) equals 
						totalDepth.  Some of the functions from Pagel, above, change the total depth of the tree, so comparing 
						parameters estimated under these models can be difficult.				
						
						
				
}
\value{
	Transformed phylogenetic tree, object of class phylo (ape format).
}
\references{Delta and lambda: Pagel, M. 1999. Inferring the historical patterns of biological evolution. Nature 401:877-884.
			OU: Butler, M.A. and A.A. King, 2004. Phylogenetic comparative analysis: a modeling approach for adaptive evolution. American Naturalist 164:683-695.
			Others: Various papers in prep., L. J. Harmon and J. T. Weir.}
\author{Luke J. Harmon}
\examples{


data(geospiza)
attach(geospiza)

deltaTree(geospiza.tree, 0.5)->g2
plot(g2)
lambdaTree(geospiza.tree, 0.5)->g3
plot(g3)
kappaTree(geospiza.tree, 0.5)->g3b
plot(g3b)
tworateTree(geospiza.tree, 0.5, 0.5)->g4
plot(g4)
linearchangeTree(geospiza.tree, 0.1)->g5
plot(g5)
exponentialchangeTree(geospiza.tree, 0.1)->g5
plot(g5)
speciationalTree(geospiza.tree)->g6
plot(g6)
rescaleTree(geospiza.tree, 100)->g7
plot(g7)


}
\keyword{ manip }
