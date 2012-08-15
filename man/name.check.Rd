\name{name.check}
\alias{name.check}
\alias{treedata}
\title{ Compares taxa in data and tree }
\description{
This function is a general tool for checking for concordance between a data file and a phylogenetic
tree.  For the data, names can be specified as the names of objects in the vector or rownames of the data 
array. The name.check function finds and lists all taxa present in data set but not 
in the tree, and vice-versa.  The treedata function returns a list containing both the tree and the data 
after pruning out any species that are not found in both. }
\usage{
name.check(phy, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{ an object of class "phylo" }
  \item{data}{ data for tips of the tree }
}
\value{
  \item{tree_not_data}{Taxa in tree but not data}
  \item{data_not_tree}{Taxa in data but not tree}
  ...
}
\author{ Luke J. Harmon }
\examples{

geo=get(data(geospiza))

name.check(geo$phy, geo$dat)

}

\keyword{ arith }
