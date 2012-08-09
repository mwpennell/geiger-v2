\name{fitDiscrete}
\alias{fitDiscrete}
\title{ Model fitting for discrete comparative data }
\description{
	Fits macroevolutionary models to phylogenetic trees
}
\usage{
fitDiscrete(phy, dat, model = c("ER", "SYM", "ARD", "meristic"),
                 transform = c("none", "EB", "lambda", "kappa","delta", "white"), 
				 bounds = list(), control = list(method = c("SANN", "L-BFGS-B"), 
				 niter = 100, FAIL = 1e+200, hessian = FALSE, hessian_P = 0.05), ...)
}
\arguments{
  \item{phy}{ a phylogenetic tree of class phylo}
  \item{dat}{ data vector for a single trait, with names matching tips in \code{phy}}
  \item{model}{ model to fit to comparative data (see \bold{Details}) }
  \item{bounds}{ range to constrain parameter estimates (see \bold{Details})  }
  \item{control}{ settings used for optimization of the model likelihood}
  \item{...}{if \code{model="meristic"}, \code{...} can dictate whether the matrix is asymmetric (\code{symmetric=FALSE})}
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

gb<-round(geospiza.data[,5])
names(gb)<-rownames(geospiza.data)
mm=matrix(sample(1:3, 9, replace=TRUE),3)
diag(mm)=NA

fitDiscrete(geospiza.tree, gb, model="ER")
fitDiscrete(geospiza.tree, gb, model="ARD")
fitDiscrete(geospiza.tree, gb, model=mm)
fitDiscrete(geospiza.tree, gb, model=mm, transform="delta")




}
\keyword{ arith }
