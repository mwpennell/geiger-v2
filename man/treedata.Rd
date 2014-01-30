\name{treedata}
\alias{treedata}
\title{ compare taxa in data and tree }
\description{
matching species found in a comparative dataset}

\usage{
treedata(phy, data, sort=FALSE, warnings=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{data}{a named vector or matrix of continuous trait values, for species in \code{phy}}
  \item{sort}{whether to sort the \code{data} based on names found in \code{phy}}
  \item{warnings}{whether to report warnings of mismatched taxa between \code{phy} and \code{data}}
}
\details{
This function is a general tool for checking for concordance between a data object and a phylogenetic
tree.  For the \code{data}, names can be specified as the names of objects in the vector or rownames of the data 
array. 
}
\value{
  The function returns a list of two elements (\code{phy} and \code{data}) that are manipulated to include only those species 
  found in both the tree and data supplied by the user.  
}
\author{ LJ Harmon }
\examples{

geo <- get(data(geospiza))

treedata(geo$phy, geo$dat, sort=TRUE, warnings=TRUE)

}

\keyword{manip}
\keyword{data}
\keyword{graphs}

