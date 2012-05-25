\name{fitDiscrete}
\alias{fitDiscrete}
\title{ Model fitting for discrete comparative data }
\description{
	Fits macroevolutionary models to phylogenetic trees
}
\usage{
fitDiscrete(phy, data, model=c("ER", "SYM", "ARD"), treeTransform=c("none", "lambda", "kappa", "delta", "linearChange", "exponentialChange", "twoRate"), data.names=NULL, plotlnl=F, qLimits=c(0.0001, 1000), pLimits=c(0.00001, 10))
}
\arguments{
  \item{phy}{ object of type phylo }
  \item{data}{ Data vector (one trait) or matrix (multiple traits) }
  \item{model}{ One of ER, SYM, or ARD; see below }
  \item{treeTransform}{ Model for transforming the tree; see below } 
  \item{data.names}{ Tip names for data vector that match tree species; ignored if data includes names}
  \item{plotlnl}{ Plot likelihood surface?  Works only with models of two parameters or fewer.}
  \item{qLimits}{Vector giving minimum and maximum values for rate parameter q; used ONLY for plotting.}
  \item{pLimits}{Vector giving minimum and maximum values for whichever tree transformation parameter you're using; used ONLY for plotting.}

}
\details{
  This function fits various likelihood models for discrete character evolution.   Likelihood is maximized 
  using the r function nlm. All of the models are 
  continuous-time Markov models of trait evolution (see Yang 2006 for a good general discussion of this 
  type of model).  The function can handle traits with any number of character states, under a range of models. 
  The character model is specified by the "model" argument:
  ER: Equal-rates (Mk) model; all transitions occur at equal rates. 
  SYM: Symmetric transitions are equal; that is, 0->1 occurs at the same rate as 1->0, which may 
  differ from the transition rate between states 1 and 2. For a 2 state character, this model is equivalent 
  to ER  
  ARD: All rates different model; each rate is a separate parameter. 
  
  The function returns a rate matrix, Q, giving the transition rates among the characters.  The diagonals of this matrix are zero; 
  this is only because the rows of the matrix sum to zero.
  
  The treeTransform argument allows you to test models where rates vary across the tree.  Options are:
 none: Rates are constant through time 
 lambda: Pagel's lambda; multiplies all internal branches of the tree by lambda, leaving tip branches as their original length.
 kappa: Pagel's kappa; raises all branch lengths to the power kappa.  As kappa approaches zero, the model becomes speciational.
 delta: Pagel's delta; raises all node depths to the power delta. If delta is less than one, evolution in concentrated early in the tree; delta > 1 concentrates evolution towards the tips.
 exponentialChange: A model where the rate of evolution changes exponentially through time.
						The difference between this option and linearChange is that rates grow or decay exponentially rather than linearly.
						The equation is r(t) = ro * exp(a * t), where ro is the inital rate and a is the rate of rate change.
 linearChange: A model where the rate of evolution changes linearly through time.
						Rates change linearly: r(t) = ro + a * t, where ro is the inital rate and a is the slope determined using endRate. 
						If endRate is greater than 1, evolution gradually speeds up; if endrate is less than 1, evolution gradually slows down. If endRate = 1, the model is a constant-
						rate model.

 twoRate: A model that effectively changes the rate of evolution at some point in time to endRate. 
 						If endRate is greater than 1, evolution speeds up, all branches towards the tips are made longer,  
						if endrate is less than 1, branches towards the tips are made shorter. If endRate = 1, the model is a constant-
						rate model.
}
\value{
 	Returns maximum likelihood value for q and selected parameters, along with the likelihood score.
	Sometimes, parameters are confounded!
}
\references{ 
Yang, Z. 2006. Computational Molecular Evolution. Oxford University Press: Oxford.
}
\author{ Luke J. Harmon and R. E. Glor }
\examples{
data(geospiza)
attach(geospiza)

gb<-as.factor(geospiza.data[,1]>4.2)
names(gb)<-rownames(geospiza.data)

fitDiscrete(geospiza.tree, gb)
fitDiscrete(geospiza.tree, gb, treeTransform="lambda")
fitDiscrete(geospiza.tree, gb, treeTransform="delta")
fitDiscrete(geospiza.tree, gb, treeTransform="kappa")

fitDiscrete(geospiza.tree, gb, treeTransform="linearChange")
fitDiscrete(geospiza.tree, gb, treeTransform="exponentialChange")
fitDiscrete(geospiza.tree, gb, treeTransform="twoRate")




}
\keyword{ arith }
