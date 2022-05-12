\name{subset.phylo}
\alias{subset.phylo}
\title{
blending information from taxonomies and trees }
\description{
working with systematic reference tables and phylogenies}
\usage{

    \method{subset}{phylo}(x, taxonomy, rank="", ncores=1, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a phylogenetic tree of class 'phylo'}
  \item{taxonomy}{a ('matrix') linkage table between tips of the phylogeny and clades represented in the tree; rownames of 'taxonomy' should be tips found in the phylogeny}
  \item{rank}{a column name in 'taxonomy' at which to resample the tree (see \bold{Examples})}
  \item{ncores}{max number of cores to use}
  \item{...}{arguments to be passed to other functions (has no effect in the present context)}
}

\author{JM Eastman}

\examples{
\dontrun{
sal <- get(data(caudata))
print(head(sal$tax))

nphy <- subset(sal$phy, sal$tax, "genus", ncores=1)
plot(nphy, type="fan", cex=0.15)
}
}

\keyword{ graphs }
\keyword{manip}
