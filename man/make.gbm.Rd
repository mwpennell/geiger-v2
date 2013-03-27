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
    The argument \code{...} controls the substitution of default settings for Markov-chain Monte Carlo sampling. Below are the settings that are controllable by the user. 
	These parameters and their default settings can also be found with an \emph{empty} call to the function (e.g., \code{make.gbm()}).
	
	\bold{measurement error (SE)}: one of the arguments necessary for running \code{\link[geiger]{rjmcmc.bm}} is \code{SE}, which is a statement about the error associated with the values 
	given in \code{dat}.  Measurement error (whose argument is \code{SE}) can be a named vector of numeric values (including \code{NA}) or a single value (including \code{NA}).  If given as 
	a vector, \code{SE} must have names that correspond to the those found for \code{dat}. If given a single value for 
	\code{SE}, the sampler will apply that value of measurement error to all tips in the tree. If \code{NA} appears for the measurement error for \emph{any} species, 
	\code{SE} becomes an additional parameter of the model: this density is consequently sampled by \code{rjmcmc.bm}. The default for \code{rjmcmc.bm} is to estimate a single \code{SE} 
	(which is applied to all species).  
	
	\bold{control settings}: default settings for each control parameter are given below.  Note that for the discrete random variables (for which \code{dlnSHIFT} and \code{dlnJUMP}) apply, 
	certain criteria must be met if the user prefers to supply a different prior density. The function \code{\link[geiger]{dcount}} is useful for building a custom prior density function 
	for discrete variables.
}
\value{
The items that can be tailored in the resulting control object are as follows:
	
				\item{\bold{method}}{default \code{direct}; this defines the implementation for computing model likelihood (currently with only the default available)}				
				\item{\bold{constrainSHIFT}}{default \code{FALSE}; if integer given, this constrains the number of local clocks in the diffusive process}
				\item{\bold{constrainJUMP}}{default \code{FALSE}; if integer given, this constrains the number of pulses in the jump process}
				\item{\bold{dlnSHIFT}}{default \code{dcount(0:(nn - 1), FUN = dpois, lambda = log(2))}, where \code{nn} is the number of branches in the tree; this controls the prior density on the number of shifts between local clocks in the tree (and applies only if the \code{type} argument to \code{\link[geiger]{rjmcmc.bm}} is \code{"rbm"} or \code{"jump-rbm"}) and if \code{constrainSHIFT} is \code{FALSE})}
				\item{\bold{dlnJUMP}}{default \code{dcount(0:nn, FUN = dlunif, min = 0, max = nn, dzero = 0.5)}, where \code{nn} is the number of branches in the tree; this controls the prior density on the number of evolutionary pulses across the tree (and applies only if the \code{type} argument to \code{\link[geiger]{rjmcmc.bm}} is \code{"jump-bm"} or \code{"jump-rbm"} and if \code{constrainJUMP} is \code{FALSE})}
				\item{\bold{dlnRATE}}{default \code{function (x) dexp(x, rate = 1/(10^3), log = TRUE)}; this defines the prior density on rate scalars}
				\item{\bold{dlnSE}}{default \code{function (x) dexp(x, rate = 1/(10^3), log = TRUE)}; this defines the prior density on measurement error}
				\item{\bold{dlnPULS}}{default \code{function (x) dexp(x, rate = 1/(10^3), log = TRUE)}; this defines the prior density on jump variance}
				\item{\bold{dlnROOT}}{default \code{function (x) dunif(x, min=-10^3, max=10^3, log=TRUE)}; this defines the prior density on root state}
				\item{\bold{rate.lim}}{default \code{list(min=0, max=Inf)}; this defines the numerical limits on the rate scalars}
				\item{\bold{se.lim}}{default \code{list(min=0, max=Inf)}; this defines the numerical limits on the measurement error}
				\item{\bold{root.lim}}{default \code{list(min=-Inf, max=Inf)}; this defines the numerical limits on the root state}
				\item{\bold{jump.lim}}{default \code{1}; determines the number of jumps permissible along each branch (1 is currently supported)}
				\item{\bold{excludeSHIFT}}{default \code{c()}; if this argument is not empty, these are branches (specified by numeric node-identifiers) which cannot be chosen for a novel local clock}
				\item{\bold{excludeJUMP}}{default \code{c()}; if this argument if not empty, these are branches (specified by numeric node-identifiers) which cannot be chosen for an evolutionary pulse}
				\item{\bold{bm.jump}}{default \code{0.5}; this defines the ratio between proposals for the diffusion and jump processes (0.5 is perfectly balanced)}
				\item{\bold{mergesplit.shift}}{default \code{0.5}; this defines the ratio between proposals that scale model complexity and those that do not}
				\item{\bold{tune.scale}}{default \code{0.65}; this defines the ratio between proposals that modify the model branchwise versus treewide}
				\item{\bold{slide.mult}}{default \code{0.25}; this defines the ratio between sliding window and multiplier proposals}
				\item{\bold{prob.dimension}}{default \code{0.65}; this defines the proportion of proposals used to modify model dimensionality}
				\item{\bold{prob.effect}}{default \code{0.3}; this defines the proportion of proposals that do not alter dimensionality}
				\item{\bold{prob.SE}}{default \code{0.03}; this defines the proportion of proposals that alter the (unknown) measurement error}
				\item{\bold{prob.root}}{default \code{0.02}; this defines the proportion of proposals that alter the root state}
				\item{\bold{prop.width}}{default \code{1}; this defines the proposal width used for multiplier and sliding-window proposals}
				\item{\bold{simple.start}}{default \code{TRUE}; this determines whether to start the MCMC sampler with minimal dimensionality}
				\item{\bold{filebase}}{default \code{"result"}; this defines a 'base' filename for the output}
		 
}

\author{JM Eastman}

\seealso{
	\code{\link[geiger]{rjmcmc.bm}}
}

\keyword{models}


