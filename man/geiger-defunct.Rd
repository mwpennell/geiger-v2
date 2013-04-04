\name{geiger-defunct}
\alias{area.between.curves}
\alias{BDsim}
\alias{birthdeath.tree}
\alias{calibrate.proposalwidth}
\alias{compare.rates}
\alias{deltaTree}
\alias{disp.calc}
\alias{dtt.full}
\alias{exponentialchangeTree}
\alias{get.simulation.matrix}
\alias{getAncStates}
\alias{ic.sigma}
\alias{intercalate.samples}
\alias{kappaTree}
\alias{lambdaTree}
\alias{linearchangeTree}
\alias{name.check}
\alias{node.leaves}
\alias{node.sons}
\alias{ouTree}
\alias{phy.anova}
\alias{phy.manova}
\alias{pool.rjmcmcsamples}
\alias{prune.extinct.taxa}
\alias{prune.random.taxa}
\alias{rate.estimate}
\alias{rescaleTree}
\alias{runMedusa}
\alias{shifts.plot}
\alias{speciationalTree}
\alias{tip.disparity}
\alias{tracer}
\alias{tworateTree}
\alias{vmat}
\title{deprecated functions in GEIGER}
\description{
  This function has been deprecated in \pkg{geiger} and \pkg{auteur} or may be available from another package. Below shows the original function and the suggested function to use in its place.
}
\details{
	\itemize{
		\item{\code{area.between.curves}}: use \code{geiger:::.area.between.curves}
		\item{\code{BDsim}}: use \code{\link[geiger]{sim.bd}}
		\item{\code{birthdeath.tree}}: use \code{\link[geiger]{sim.bdtree}}
		\item{\code{calibrate.proposalwidth}}: use \code{\link[geiger]{calibrate.rjmcmc}}
		\item{\code{deltaTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{disp.calc}}: use \code{\link[geiger]{disparity}}
		\item{\code{dtt.full}}: use \code{\link[geiger]{dtt}}
		\item{\code{exponentialchangeTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{get.simulation.matrix}}: use \code{geiger:::.get.simulation.matrix}
		\item{\code{getAncStates}}: use \code{\link[phytools]{fastAnc}}
		\item{\code{ic.sigma}}: use \code{\link[geiger]{vcv.phylo}}
		\item{\code{intercalate.samples}}: use \code{geiger:::.intercalate.rjmcmc}
		\item{\code{kappaTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{lambdaTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{linearchangeTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{name.check}}: use \code{geiger:::.treedata}
		\item{\code{node.leaves}}: use \code{\link[geiger]{tips}}
		\item{\code{node.sons}}: use \code{geiger:::.get.desc.of.node}
		\item{\code{ouTree}}: use \code{\link[geiger]{transform.phylo}}
        \item{\code{phy.anova}}: use \code{\link[geiger]{aov.phylo}}
        \item{\code{phy.manova}}: use \code{\link[geiger]{aov.phylo}}
		\item{\code{pool.rjmcmcsamples}}: use \code{\link[geiger]{load}}
		\item{\code{prune.extinct.taxa}}: use \code{\link[geiger]{drop.extinct}}
		\item{\code{prune.random.taxa}}: use \code{\link[geiger]{drop.random}}
		\item{\code{rate.estimate}}: use \code{\link[geiger]{bd.ms}} or \code{\link[geiger]{bd.km}}
		\item{\code{rescaleTree}}: use \code{\link[geiger]{transform.phylo}}
        \item{\code{runMedusa}}: use \code{\link[geiger]{medusa}}
		\item{\code{shifts.plot}}: use \code{plot}
		\item{\code{speciationalTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{tip.disparity}}: use \code{\link[geiger]{disparity}}
		\item{\code{tracer}}: use of \code{\link[coda]{plot.mcmc}} and other functions in \pkg{coda} recommended
		\item{\code{tworateTree}}: use \code{\link[geiger]{transform.phylo}}
		\item{\code{vmat}}: use \code{geiger:::.vmat}
	}
}

