\name{vcv.phylo}
\alias{vcv.phylo}
\title{evolutionary or phylogenetic VCV matrix}
\description{
	estimating the evolutionary or phylogenetic variance-covariance matrix}
\usage{
\method{vcv}{phylo}(phy, data=NULL, \dots)
}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{data}{a named vector or matrix of continuous trait values, associated with species in \code{phy}}
  \item{...}{arguments to be passed to other functions (has no effect in the present context)}
}
\details{
	If given \code{data} for \code{n} quantitative variables, this function returns the estimated evolutionary variance-covariance matrix of the variables under a multivariate Brownian motion model.  Note that other evolutionary models may be possible if the tree is first transformed (see \code{\link[geiger]{transform.phylo}} and \bold{Examples}).  If you have \code{n} characters in your analysis, this will be an \code{n}x\code{n} matrix. Diagonal elements represent rate estimates for individual characters, while off-diagonal elements represent the estimated covariance between two characters.
	
	If the \code{data} object is \code{NULL}, a phylogenetic variance-covariance matrix is returned (as in \code{\link[ape]{vcv.phylo}} of \pkg{ape}). The result of this function is the \code{N}-by-\code{N} matrix of expected trait variances and trait covariances among species under Brownian-motion evolution (for a tree with \code{N} species).  

}
\references{ Revell, L. J., L. J. Harmon, R. B. Langerhans, and J. J. Kolbe. 2007. A phylogenetic 
             approach to determining the importance of constraint on phenotypic evolution in the 
             neotropical lizard, \emph{Anolis cristatellus}. \emph{Evolutionary Ecology Research} 9: 261-282. }
\author{ LJ Harmon }
\seealso{
	\code{\link[ape]{vcv.phylo}} in \pkg{ape}
}
\examples{

geo=get(data(geospiza))

## EVOLUTIONARY VCV
vcv(geo$phy, geo$dat)

## EVOLUTIONARY VCV -- assuming speciational model
kphy=transform(geo$phy, "kappa", 0)
vcv(kphy, geo$dat)

## PHYLOGENETIC VCV
vcv(geo$phy)



}

\keyword{ arith }
