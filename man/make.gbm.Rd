\name{make.gbm}
\alias{make.gbm}
\title{
tailor reversible-jump Markov chain Monte Carlo sampling}
\description{
controlling reversible-jump Markov chain Monte Carlo sampling}
\usage{
make.gbm(phy, dat, SE=NA, type = c("bm", "rbm", "jump-bm", "jump-rbm"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{dat}{a named vector of continuous trait values, associated with each species in \code{phy}}
  \item{SE}{a named vector of standard errors for each trait value; applied to all trait values if given a single value}
  \item{type}{the class of model to use (see \code{\link[geiger]{rjmcmc.bm}})}
  \item{...}{arguments passed internally to control other settings (see \bold{Details})}

}
\details{
	Currently, the only \code{type} possible is for the generalized model for continuous traits (i.e., \code{type="rbm"}).  
	
	The argument \code{...} controls the substitution of default settings for Markov-chain Monte Carlo sampling. Below are the settings that are controllable by the user. 
	These parameters and their default settings can be found with an \emph{empty} call to the function (e.g., \code{make.gbm()}).
	
	\bold{measurement error (SE)}: one of the arguments necessary for running \code{\link[geiger]{rjmcmc.bm}} is \code{SE}, which is a statement about the error associated with the values 
	given in \code{dat}.  Measurement error (whose argument is \code{SE}) can be a named vector of numeric values (including \code{NA}) or a single value (including \code{NA}).  If given as 
	a vector, \code{SE} must have names that correspond to the those found for \code{dat}. If given a single value for 
	\code{SE}, the sampler will apply that value of measurement error to all tips in the tree. If \code{NA} appears for the measurement error for \emph{any} species, 
	\code{SE} becomes an additional parameter of the model: this density is consequently sampled by \code{rjmcmc.bm}. The default for \code{rjmcmc.bm} is to estimate a single \code{SE} 
	(which is applied to all species).  
	
	\bold{control settings}: default settings for each control parameter are given below.  Note that for the discrete random variables (for which \code{dlnSHIFT} and \code{dlnJUMP}) apply, 
	certain criteria must be met if the user prefers to supply a different prior density. The function \code{\link[geiger]{dcount}} is strongly recommended for building a custom prior density function 
	for discrete variables.
	
			\itemize{
				\item{\emph{method} -- }{\code{direct}; the implementation for computing model likelihood}				
				\item{\emph{constrainSHIFT} -- }{\code{FALSE}; if integer given, this constrains the number of local clocks in the diffusive process}
				\item{\emph{constrainJUMP} -- }{\code{FALSE}; if integer given, this constrains the number of pulses in the jump process}
				\item{\emph{dlnSHIFT} -- }{\code{dcount(0:(nn - 1), FUN = dpois, lambda = log(2))}, where \code{nn} is the number of branches in the tree}
				\item{\emph{dlnJUMP} -- }{\code{dcount(0:nn, FUN = dlunif, min = 0, max = nn, dzero = 0.5)}, where \code{nn} is the number of branches in the tree}
				\item{\emph{dlnRATE} -- }{\code{function (x) dexp(x, rate = 1/(10^3), log = TRUE)}; prior density on rate scalars}
				\item{\emph{dlnSE} -- }{\code{function (x) dexp(x, rate = 1/(10^3), log = TRUE)}; prior density on measurement error}
				\item{\emph{dlnPULS} -- }{\code{function (x) dexp(x, rate = 1/(10^3), log = TRUE)}; prior density on jump variance}
				\item{\emph{dlnROOT} -- }{\code{function (x) dunif(x, min=-10^3, max=10^3, log=TRUE)}; prior density on root state}
				\item{\emph{rate.lim} -- }{\code{list(min=0, max=Inf)}; limits on the rate scalars}
				\item{\emph{se.lim} -- }{\code{list(min=0, max=Inf)}; limits on the measurement error}
				\item{\emph{root.lim} -- }{\code{list(min=-Inf, max=Inf)}; limits on the root state}
				\item{\emph{jump.lim} -- }{\code{1}; determines the number of jumps permissible along each branch (1 is currently supported)}
				\item{\emph{excludeSHIFT} -- }{\code{c()}; if not empty, these are branches (specified by numeric node-identifiers) which cannot be chosen for a novel local clock}
				\item{\emph{excludeJUMP} -- }{\code{c()}; if not empty, these are branches (specified by numeric node-identifiers) which cannot be chosen for an evolutionary pulse}
				\item{\emph{bm.jump} -- }{\code{0.5}; the ratio between proposals for the diffusion and jump processes (0.5 is perfectly balanced)}
				\item{\emph{mergesplit.shift} -- }{\code{0.5}; the ratio between proposals that scale model complexity and those that do not}
				\item{\emph{tune.scale} -- }{\code{0.65}; the ratio between proposals that modify the model branchwise versus treewide}
				\item{\emph{slide.mult} -- }{\code{0.25}; the ratio between sliding window and multiplier proposals}
				\item{\emph{prob.dimension} -- }{\code{0.65}; the proportion of proposals used to modify model dimensionality}
				\item{\emph{prob.effect} -- }{\code{0.3}; the proportion of proposals that do not alter dimensionality}
				\item{\emph{prob.SE} -- }{\code{0.03}; the proportion of proposals that alter the (unknown) measurement error}
				\item{\emph{prob.root} -- }{\code{0.02}; the proportion of proposals that alter the root state}
				\item{\emph{prop.width} -- }{\code{1}; the proposal width used for multiplier and sliding-window proposals}
				\item{\emph{simple.start} -- }{\code{TRUE}; whether to start the MCMC sampler with minimal dimensionality}
				\item{\emph{filebase} -- }{\code{result}; a character string specifying a 'base' for the output}
		}
		 
}

\author{JM Eastman}

\seealso{
	\code{\link[geiger]{rjmcmc.bm}})
}


