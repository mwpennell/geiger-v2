\name{auteur-defunct}
\alias{calibrate.proposalwidth}
\alias{pool.rjmcmcsamples}
\alias{intercalate.samples}
\alias{tracer}
\alias{shifts.plot}
\alias{compare.rates}
\alias{vmat}
\alias{disp.calc}
\alias{dtt.full}
\alias{name.check}
\alias{node.sons}
\alias{prune.extinct.taxa}
\alias{prune.random.taxa}
\alias{get.simulation.matrix}
\alias{tip.disparity}
\alias{deltaTree}
\alias{lambdaTree}
\alias{kappaTree}
\alias{ouTree}
\alias{tworateTree}
\alias{linearchangeTree}
\alias{exponentialchangeTree}
\alias{speciationalTree}
\alias{rescaleTree}
\alias{getAncStates}
\alias{area.between.curves}
\alias{ic.sigma}
\alias{rate.estimate}
\alias{node.leaves}
\alias{BDsim}
\alias{birthdeath.tree}
\title{deprecated functions in GEIGER}
\description{
  This function has been deprecated in \pkg{geiger} or moved to another package.
}
\details{
	\itemize{
		\item{\code{calibrate.proposalwidth}}: use \code{\link[geiger]{calibrate.rjmcmc}}
		\item{\code{pool.rjmcmcsamples}}: use \code{\link[geiger]{load}}
		\item{\code{tracer}}: use of \code{\link[coda]{plot.mcmc}} and other functions in \pkg{coda} recommended
		\item{\code{intercalate.samples}}: use \code{geiger:::.intercalate.rjmcmc}
		\item{\code{shifts.plot}}: use \code{plot}
		\item{\code{vmat}}: use \code{geiger:::.vmat}
		\item{\code{disp.calc}}: use \code{\link[geiger]{disparity}}
		\item{\code{dtt.full}}: use \code{\link[geiger]{dtt}}
		\item{\code{name.check}}: use \code{geiger:::.treedata}
		\item{\code{node.sons}}: use \code{geiger:::.get.desc.of.node}
		\item{\code{prune.random.taxa}}: use \code{\link[geiger]{drop.random}}
		\item{\code{prune.extinct.taxa}}: use \code{\link[geiger]{drop.extinct}}
		\item{\code{get.simulation.matrix}}: use \code{geiger:::.get.simulation.matrix}
		\item{\code{tip.disparity}}: use \code{\link[geiger]{disparity}}
		\item{\code{deltaTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{lambdaTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{kappaTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{ouTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{tworateTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{linearchangeTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{exponentialchangeTree}}: \code{\link[geiger]{transform.phylo}}
		\item{\code{speciationalTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{rescaleTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{getAncStates}}: use \code{geiger:::.ace}
		\item{\code{area.between.curves}}: use \code{geiger:::.area.between.curves}
		\item{\code{ic.sigma}}: use \code{\link[geiger]{vcv.phylo}}
		\item{\code{rate.estimate}}: use \code{\link[geiger]{bd.ms}} or \code{\link[geiger]{bd.km}}
		\item{\code{node.leaves}}: use \code{\link[geiger]{tips}}
		\item{\code{BDsim}}: use \code{\link[geiger]{sim.bd}}
		\item{\code{birthdeath.tree}}: use \code{\link[geiger]{sim.bdtree}}
	}
}

\keyword{internal}
