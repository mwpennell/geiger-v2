\name{cache}
\alias{cache}
\title{
Run and parameter control for reversible-jump Markov chain Monte Carlo sampling}
\description{
Run and parameter control for reversible-jump Markov chain Monte Carlo sampling}
\usage{
cache(phy, dat, type = "bm", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{dat}{a named vector of continuous trait values, associated with each species in \code{phy}}
  \item{type}{currently the only supported option is \code{bm}}
  \item{...}{arguments passed internally to control other settings (see \bold{Details})}

}
\details{
	Currently, the only \code{type} possible is for the generalized model for continuous traits (i.e., \code{type="bm"}).  
	
	The argument \code{...} controls the substitution of default settings for Markov-chain Monte Carlo sampling. Below are the settings that are controllable by the user.
	
	\bold{measurement error (SE)}: one of the arguments necessary for running \code{\link[geiger]{rjmcmc.bm}} is \code{SE}, which is a statement about the error associated with the values 
	given in \code{dat}.  Measurement error (whose argument is \code{SE}) can be a named vector of numeric values (including \code{NA}) or a single value (including \code{NA}).  If given as 
	a vector, \code{SE} must have names that correspond to the those found for \code{dat}. If given a single value for 
	\code{SE}, the sampler will apply that value of measurement error to all tips in the tree. If \code{NA} appears for the measurement error for \emph{any} species, 
	\code{SE} becomes an additional parameter of the model: this density is consequently sampled by \code{rjmcmc.bm}. The default for \code{rjmcmc.bm} is to estimate a single \code{SE} 
	(which is applied to all species).  
	
	\bold{control settings}: default settings for each control parameter are given below.  Note that for the discrete random variables (for which \code{dlnSHIFT} and \code{dlnJUMP}) apply, 
	certain criteria must be met if the user prefers to supply a different prior density. The function \code{\link[geiger]{dcount}} is strongly recommended for building a custom prior density function 
	for discrete variables.)
	
			\itemize{
				\item{method -- }{\code{direct}; the implementation for computing model likelihood}				
				\item{constrainSHIFT -- }{\code{FALSE}; if integer given, this constrains the number of local clocks in the diffusive process}
				\item{constrainJUMP -- }{\code{FALSE}; if integer given, this constrains the number of pulses in the jump process}
				\item{dlnSHIFT -- }{\code{dcount(0:(nn - 1), FUN = dpois, lambda = log(2))}, where \code{nn} is the number of branches in the tree}
				\item{dlnJUMP -- }{\code{dcount(0:nn, FUN = dlunif, min = 0, max = nn, dzero = 0.5)}, where \code{nn} is the number of branches in the tree}
				\item{dlnRATE -- }{\code{function (x) dexp(x, rate = 1/(10^3), log = TRUE)}; prior density on rate scalars}
				\item{dlnSE -- }{\code{function (x) dexp(x, rate = 1/(10^3), log = TRUE)}; prior density on measurement error}
				\item{dlnPULS -- }{\code{function (x) dexp(x, rate = 1/(10^3), log = TRUE)}; prior density on jump variance}
				\item{dlnROOT -- }{\code{function (x) dunif(x, min=-10^3, max=10^3, log=TRUE)}; prior density on root state}
				\item{rate.lim -- }{\code{list(min=0, max=Inf)}; limits on the rate scalars}
				\item{se.lim -- }{\code{list(min=0, max=Inf)}; limits on the (unknown) measurement error}
				\item{root.lim -- }{\code{list(min=-Inf, max=Inf)}; limits on the root state}
				\item{jump.lim -- }{\code{1}; determines the number of jumps permissible along each branch (1 is currently supported)}
				\item{excludeSHIFT -- }{\code{c()}; if not empty, these are branches (specified by numeric ID) which cannot be chosen for a novel local clock}
				\item{excludeJUMP -- }{\code{c()}; if not empty, these are branches (specified by numeric ID) which cannot be chosen for an evolutionary pulse}
				\item{bm.jump -- }{\code{0.5}; the ratio between proposals for the diffusion and jump processes (0.5 is perfectly balanced)}
				\item{mergesplit.shift -- }{\code{0.5}; the ratio between proposals that scale model complexity and those that do not}
				\item{tune.scale -- }{\code{0.65}; the ratio between proposals that modify the model branchwise versus treewide}
				\item{slide.mult -- }{\code{0.25}; the ratio between sliding window and multiplier proposals}
				\item{prob.dimension -- }{\code{0.65}; the proportion of proposals used to modify model dimensionality}
				\item{prob.effect -- }{\code{0.3}; the proportion of proposals that do not alter dimensionality}
				\item{prob.SE -- }{\code{0.03}; the proportion of proposals that alter the (unknown) measurement error}
				\item{prob.root -- }{\code{0.02}; the proportion of proposals that alter the root state}
				\item{prop.width -- }{\code{1}; the proposal width used for multiplier and sliding-window proposals}
				\item{simple.start -- }{\code{TRUE}; whether to start the MCMC sampler with minimal dimensionality}
				\item{hashtips -- }{\code{c()}; a character vector of names (e.g., \code{phy$tip.label}) from which unique branch identifiers will be constructed}
				\item{filebase -- }{\code{result}; a character string specifying a 'base' for the output}
		}
		 
}

\author{JM Eastman}


