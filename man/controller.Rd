\name{controller}
\alias{controller}
\title{
Run and parameter control for reversible-jump Markov chain Monte Carlo sampling}
\description{
Run and parameter control for reversible-jump Markov chain Monte Carlo sampling}
\usage{
controller(phy, dat, SE = NULL, type = "bm", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{dat}{a named vector of continuous trait values, associated with each species in \code{phy}}
  \item{SE}{a named vector of standard errors for each trait value; applied to all trait values if given a single value}
  \item{type}{currently the only supported option is \code{bm}}
  \item{...}{arguments passed internally to control other settings (see \bold{Details})}

}
\details{
	Default settings for each control parameter are given below.  Note that for the discrete random variables (for which \code{dlnSHIFT} and \code{dlnJUMP}) apply, certain criteria 
	must be met if the user prefers to supply a different prior density (and the function \code{\link[auteur]{dcount}} is strongly recommended for building a prior density function 
	for discrete variables.)
	 
	\itemize{
		\item{method }{\code{direct}; the implementation for computing model likelihood}
		\item{constrainSHIFT }{\code{FALSE}; if integer given, this constrains the number of local clocks in the diffusive process}
		\item{constrainJUMP }{\code{FALSE}; if integer given, this constrains the number of pulses in the jump process}
		\item{dlnSHIFT }{\code{dcount(0:(nn - 1), FUN = dpois, lambda = log(2))}, where \code{nn} is the number of branches in the tree}
		\item{dlnJUMP }{\code{dcount(0:nn, FUN = dlunif, min = 0, max = nn, dzero = 0.5)}, where \code{nn} is the number of branches in the tree}
		\item{dlnRATE }{\code{function (x) dexp(x, rate = 1/(10^3), log = TRUE)}}
		\item{dlnPULS }{\code{function (x) dexp(x, rate = 1/(10^3), log = TRUE)}}
		\item{dlnROOT }{\code{function (x) dunif(x, min=-10^3, max=10^3, log=TRUE)}}
		\item{rate.lim }{\code{list(min=0, max=Inf)}}
		\item{root.lim }{\code{list(min=-Inf, max=Inf)}}
		\item{jump.lim }{\code{1}; determines the number of jumps permissible along each branch (only 1 is currently supported)}
		\item{excludeSHIFT }{\code{c(); branches (specified by numeric ID) which cannot be chosen for a novel local clock}}
		\item{excludeJUMP }{\code{c(); branches (specified by numeric ID) which cannot be chosen for an evolutionary pulse}}
		\item{bm.jump }{\code{0.5}; the ratio between proposals for the diffusion and jump processes (0.5 is perfectly balanced)}
		\item{mergesplit.shift }{\code{0.5}; the ratio between proposals that scale model complexity and those that do not}
		\item{tune.scale }{\code{0.65}; the ratio between proposals that modify the model branchwise versus treewide}
		\item{slide.mult }{\code{0.25}; the ratio between sliding window and multiplier proposals}
		\item{prob.dimension }{\code{0.65}; the proportion of proposals used to modify model dimensionality}
		\item{prob.effect }{\code{0.3}; the proportion of proposals that do not alter dimensionality}
		\item{prob.root }{\code{0.05}; the proportion of proposals that alter the root state}
		\item{prop.width }{\code{1}; the proposal width used for multiplier and sliding-window proposals}
		\item{simple.start }{\code{TRUE}; whether to start the MCMC sampler with minimal dimensionality}
		\item{hashtips }{\code{c()}; a character vector of names (e.g., \code{phy$tip.label}) from which unique branch identifiers will be constructed}
		\item{filebase }{\code{result}; a character string specifying a 'base' for the output}
	}
}

\author{JM Eastman}


