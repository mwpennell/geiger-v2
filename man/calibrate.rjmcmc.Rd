\name{calibrate.rjmcmc}
\alias{calibrate.rjmcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
initialize proposal width
}
\description{
estimating a reasonable proposal width to initiate sampling for Markov sampling
}
\usage{
calibrate.rjmcmc(phy, dat, nstep = 10000, widths = 2^(-3:3), type=c("bm", 
    "rbm", "jump-bm", "jump-rbm"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{dat}{a named vector of continuous trait values, associated with each species in \code{phy}}
  \item{nstep}{number of proposal steps over which to assess proposal widths}
  \item{widths}{if unspecified, a series of proposal widths from 1/8 to 8 will be considered}
  \item{type}{a model type available in \code{\link[geiger]{rjmcmc.bm}}}
  \item{...}{arguments to be passed to \code{\link[geiger]{make.gbm}} and \code{\link[geiger]{rjmcmc.bm}}}
}
\details{
This function may be useful for constraining subsequent runs after an adequate proposal width has been approximated. MCMC samples from this calibration are not stored and do not become available to the user. This function is solely used to give the user a sense of acceptance rates that can be expected for different proposal widths. The narrower the width, the more easily the chain can become stuck. With a wider width, the chain will more quickly explore a broader parameter space, yet acceptance rates may become unacceptably low.}

\author{JM Eastman}

\seealso{
\code{\link[geiger]{rjmcmc.bm}}
}

\examples{
\donttest{
n=40
phy=rcoal(n=n)
dat=rTraitCont(phy=phy, model="BM", sigma=sqrt(0.1))
r=paste(sample(letters,9,replace=TRUE),collapse="")

## calibrate proposal width
calibrate.rjmcmc(phy=phy, dat=dat, nstep=500, widths=2^(-3:0), type="rbm")	
}
}
\keyword{data}

