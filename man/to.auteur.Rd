\name{to.auteur}
\alias{to.auteur}
\alias{to.coda}
\title{
conversion of MCMC samples between auteur and coda}
\description{
converting MCMC samples between auteur and coda}
\usage{
to.auteur(obj, phy = NULL, ...)
to.coda(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{for \code{to.auteur} -- an object of class \code{codaMCMCMC}, \code{mcmc.list}, or an object exported to an \code{rda} file by \pkg{auteur};
  for \code{to.coda} -- an object of class \code{auteurMCMCMC} or a list of objects individually of class \code{auteurMCMC}}
  \item{phy}{a phylogenetic tree of class 'phylo' against which to compile results; if NULL, the tree stored within the \code{rda} file is used}
  \item{...}{arguments (\code{burnin} and \code{thin}) to be passed to \code{\link[geiger]{load.rjmcmc}}}

}
\details{
A \code{coda} format of run(s) is recommended for diagnostic purposes; for summarization, \code{auteur} formats are advised. For single chains, the format adopted by both \pkg{auteur} and \pkg{coda} is identical (an object of class \code{mcmc}).  
For a series of combined runs, formats differ between the \pkg{auteur} and \pkg{coda} packages: \pkg{auteur} requires an intercalated (single) matrix of values, whereas functions within \pkg{coda} 
expect the values to be concatenated into a list (of class \code{mcmc.list}).  The function \code{to.coda} is used solely for pooling multiple runs into a format compatible with the \pkg{coda} package.    
}
\value{
For \code{to.auteur}, an object of class \code{auteurMCMCMC} (given multiple runs) or \code{auteurMCMC} (given a single run) is returned; 
for \code{to.coda}, an object of class \code{codaMCMCMC} is returned.
}

\seealso{\code{\link[geiger]{load.rjmcmc}}}

\author{JM Eastman}
\keyword{data}
\keyword{manip}

