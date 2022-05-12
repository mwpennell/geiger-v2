\name{rjmcmc.bm}
\alias{rjmcmc.bm}
\title{Bayesian sampling of shifts in trait evolution: relaxed Brownian motion}
\description{Implements reversible-jump Markov chain Monte Carlo sampling for trait evolutionary models}
\usage{
rjmcmc.bm(phy, dat, SE=NA, ngen = 50000, samp = 100, 
    type = c("jump-rbm", "rbm", "jump-bm", "bm"), ...)
}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{dat}{a named vector of continuous trait values, associated with each species in \code{phy}}
  \item{SE}{a named vector of standard errors for each trait value; applied to all trait values if given a single value}
  \item{ngen}{number of sampling generations}
  \item{samp}{frequency with which Markov samples are retained (e.g., \code{samp=10} retains every tenth sample in the chain)}
  \item{type}{the class of model to use (see \bold{Details})}
  \item{...}{arguments passed to \code{\link[geiger]{make.gbm}}}
}

\details{
Implemented is an MCMC sampler for a general model of Brownian motion, which in the full model (\code{type="jump-rbm"}) allows relaxed local clocks and also a point process of pulses 
in evolutionary rate along individual branches. Restricted models include global-rate Brownian motion (\code{type="bm"}), relaxed-rates Brownian motion (\code{type="rbm"}), and models including
jumps but a single rate of diffusion across the tree (\code{type="jump-bm"}). 

Where applicable, posterior estimates of \code{shifts} between local rates, estimates of the \code{rates} themselves, and inferred \code{jumps} (or pulses) are provided as output. 
Estimates are stored as an MCMC-generations-by-branches matrix (see \bold{Examples}), and branches are uniquely labeled by a cryptographic function to ensure comparability amongst trees 
differing in topology (see \code{\link[digest]{digest}}). 

Note that default settings (as the user assumes if nothing is specified in \code{...}) provide absolutely no guarantee of the chain achieving convergence.  The user is 
emphatically encouraged to supply informed arguments for what are the most critical aspects of this MCMC sampler (see \code{\link[geiger]{make.gbm}} for more information on 
permissible modifications to the MCMC sampler). Finding reasonable run parameters will likely require much trial and error.  Run diagnosis and inspection of chain mixing is 
facilitated by the R-package \pkg{coda} or by the Java application, Tracer (\url{http://tree.bio.ed.ac.uk/software/tracer/}).

In the \bold{Examples} below, do \emph{not} expect such short chains to reach stationarity!

}
\value{
After a run has completed, acceptance rates for the primary proposal mechanisms are printed to the console, along with 
settings of control parameters for the run (see \code{\link[geiger]{make.gbm}}).  

Posterior results are written to several files within a base directory, the contents of which are as follows:
  \item{\bold{log} }{is a logfile including the following for each Markov chain: the generations at which samples were retained (\code{state}), the \code{min}, \code{max}, and \code{median} rate of the 
  diffusion process across the tree, the number of evolutionary pulses (\code{jumps}) along single branches, the variance associated with the jump process (\code{jumpvar}), the \code{root} state, 
  and the likelihood (\code{lnL}) and prior (\code{lnLp}) of the model at sampled generations.}
  
  \item{\bold{rda} }{is a compressed R object, which stores branchwise estimates of the jump and diffusion processes. In order to be interpretable, the \code{rda} file should be processed by the function 
  \code{\link[geiger]{load.rjmcmc}}.  The package \pkg{coda} can be used within R for run diagnostics on the processed output (see, e.g., \code{\link[coda]{heidel.diag}} and \code{\link[coda]{autocorr}}.}
}

\author{JM Eastman, LJ Harmon, AL Hipp, and JC Uyeda}

\references{ Eastman JM, ME Alfaro, P Joyce, AL Hipp, and LJ Harmon. 2011. A novel comparative method for identifying shifts in the rate of character evolution on trees. \emph{Evolution} 65:3578-3589.
}

\seealso{
\code{\link[geiger]{load.rjmcmc}}
}

\examples{
\donttest{
## GENERATE DATA: jump-diffusion
phy <- ladderize(sim.bdtree(n=200), right=FALSE)
r <- paste(sample(letters,9,replace=TRUE),collapse="")
defpar <- par(no.readonly=TRUE)


tmp <- ex.jumpsimulator(phy, jumps=10)
dat <- tmp$dat
hist <- tmp$hist

ex.traitgram(phy, hist, alpha=0) # plot history of trait change

## RUN ANALYSIS


## coda package is not a dependency of geiger
## but is very useful for evaluating mcmc runs
## library(coda)


rjmcmc.bm(phy,dat, prop.width=1.5, ngen=20000, samp=500, filebase=r,
          simple.start=TRUE, type="jump-bm")
outdir <- paste("jump-BM", r, sep=".")
ps <- load.rjmcmc(outdir)

dev.new()
plot(x=ps, par="jumps", burnin=0.25, legend=FALSE, show.tip=FALSE, type="fan", edge.width=2)
mm=match(phy$edge[,2],hist$descendant)
hist=hist[mm,]
edgelabels.auteur(text=NULL, pch=21, cex=hist$cex, bg=NA, col=ifelse(hist$cex>0, 1, NA), lty=2)
title("red (estimated); black (true jump size)", line=-5)
par(defpar)

dev.new()
## from the coda package
coda::autocorr.plot(ps$log, ask=dev.interactive())
plot(ps$log, ask=dev.interactive())

## GENERATE DATA: multi-rate diffusion
scl <- ex.ratesimulator(phy, min=12, show.tip=FALSE)
dat <- rTraitCont(scl)

## RUN ANALYSIS
rjmcmc.bm(phy, dat, prop.width=1.5, ngen=20000, samp=500, filebase=r, simple.start=TRUE, type="rbm")
outdir <- paste("relaxedBM", r, sep=".")
ps <- load.rjmcmc(outdir)
dev.new()
plot(x=ps, par="shifts", burnin=0.25, legend=TRUE, show.tip=FALSE, edge.width=2)

if(!interactive()) { ## clean up
    dirs <- dir(pattern="^(jump-BM|relaxedBM)")
    unlink(dirs, recursive=TRUE)
}
}
}
\keyword{models}

