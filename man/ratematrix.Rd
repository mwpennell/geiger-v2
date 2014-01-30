\name{ratematrix}
\alias{ratematrix}
\title{evolutionary VCV matrix}
\description{
	estimating the evolutionary or phylogenetic variance-covariance matrix}
\usage{
ratematrix(phy, dat)
}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{dat}{a named vector or matrix of continuous trait values, associated with species in \code{phy}}
}
\details{
	If given \code{dat} for \code{n} quantitative variables, this function returns the estimated evolutionary variance-covariance matrix of the variables under a multivariate Brownian motion model.  Note that other evolutionary models may be possible if the tree is first transformed (see \code{\link[geiger]{rescale.phylo}} and \bold{Examples}).  If you have \code{n} characters in your analysis, this will be an \code{n}x\code{n} matrix. Diagonal elements represent rate estimates for individual characters, while off-diagonal elements represent the estimated covariance between two characters.
}
\references{ Revell, L. J., L. J. Harmon, R. B. Langerhans, and J. J. Kolbe. 2007. A phylogenetic 
             approach to determining the importance of constraint on phenotypic evolution in the 
             neotropical lizard, \emph{Anolis cristatellus}. \emph{Evolutionary Ecology Research} 9: 261-282. }
\author{ LJ Harmon }
\examples{

geo <- get(data(geospiza))

## EVOLUTIONARY VCV
ratematrix(geo$phy, geo$dat)

## EVOLUTIONARY VCV -- assuming speciational model
kphy <- rescale(geo$phy, "kappa", 0)
ratematrix(kphy, geo$dat)
geo <- get(data(geospiza))

## EVOLUTIONARY VCV
ratematrix(geo$phy, geo$dat)

## EVOLUTIONARY VCV -- assuming speciational model
kphy <- rescale(geo$phy, "kappa", 0)
ratematrix(kphy, geo$dat)

}

\keyword{ arith }
