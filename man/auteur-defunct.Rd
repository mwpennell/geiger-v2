\name{auteur-defunct}
\alias{calibrate.proposalwidth}
\alias{pool.rjmcmcsamples}
\alias{intercalate.samples}
\alias{tracer}
\alias{shifts.plot}
\alias{compare.rates}
\alias{vmat}
\title{deprecated functions in AUTEUR}
\description{
  This function has been deprecated in \pkg{auteur} or moved to another package.
}
\details{
	\itemize{
		\item{\code{calibrate.proposalwidth}}: use \code{\link[geiger]{calibrate.rjmcmc}}
		\item{\code{pool.rjmcmcsamples}}: use \code{\link[geiger]{load.rjmcmc}}
		\item{\code{tracer}}: use of \code{\link[coda]{plot.mcmc}} and other functions in \pkg{coda} recommended
		\item{\code{intercalate.samples}}: use \code{geiger:::.intercalate.rjmcmc}
		\item{\code{shifts.plot}}: use \code{plot}
		\item{code{vmat}}: use \code{geiger:::.vmat}
	}
}

\keyword{internal}
