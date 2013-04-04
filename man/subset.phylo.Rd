\name{subset.phylo}
\alias{subset.phylo}
\title{
blending information from taxonomies and trees }
\description{
working with systematic reference tables and phylogenies}
\usage{

    \method{subset}{phylo}(x, taxonomy, rank="", ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a phylogenetic tree of class 'phylo'}
  \item{taxonomy}{a ('matrix') linkage table between tips of the phylogeny and clades represented in the tree; rownames of 'taxonomy' should be tips found in the phylogeny}
  \item{rank}{a column name in 'taxonomy' at which to resample the tree (see \bold{Examples})}
  \item{...}{arguments to be passed to other functions (has no effect in the present context)}
}

\author{JM Eastman}

\examples{
sal=get(data(caudata))
print(head(sal$tax))

nphy=subset(sal$phy, sal$tax, "genus")
plot(nphy, type="fan", cex=0.15)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{manip}

