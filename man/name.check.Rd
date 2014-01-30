\name{name.check}
\alias{name.check}
\title{ Compares taxa in data and tree }
\description{
This function is a general tool for checking for concordance between a data file and a phylogenetic
tree.  For the data, names can be specified as the names of objects in the vector, rownames of the data 
array or as 'data.names'. The name.check function finds and lists all taxa present in data set but not 
in the tree, and vice-versa.  The treedata function returns a list containing both the tree and the data 
after pruning out any species that are not found in both. }
\usage{
name.check(phy, data, data.names=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{ an object of class "phylo" }
  \item{data}{ data for tips of the tree }
  \item{data.names}{ names of the tips in the order of the data; if this is not given, names 
  					will be taken from the names or rownames of the object data}
}
\value{
  \item{Tree.not.data}{Taxa in tree but not data}
  \item{Data.not.tree}{Taxa in data but not tree}
  ...
}
\author{ Luke J. Harmon }
\examples{

data(geospiza)

tmp <- name.check(geospiza$phy, geospiza$dat)
tmp

## then match data to tree
newphy <- drop.tip(geospiza$phy, tip=tmp$tree_not_data)

## name check should now say "OK"
name.check(newphy, geospiza$dat)

## this can all be done in one step using treedata
td <- treedata(geospiza$phy, geospiza$dat)
td

all(td$phy$tip.label == newphy$tip.label)

}

\keyword{ arith }
