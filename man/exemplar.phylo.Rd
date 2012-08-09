\name{exemplar.phylo}
\alias{exemplar.phylo}
\title{
Blending information from taxonomies and trees }
\description{
Working with systematic reference tables and phylogenies}
\usage{
exemplar.phylo(phy, strict.vals=NULL, taxonomy=NULL, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo' ('multiPhylo' in \code{phylo.clades})}
  \item{strict.vals}{if non-NULL, a character vector specifying the columns of \code{taxonomy} to be used in relabeling}
  \item{taxonomy}{a ('matrix') linkage table between tips of the phylogeny and clades represented in the tree; rownames of 'taxonomy' should be tips found in the phylogeny}
  \item{...}{additional arguments to be passed to \code{\link[geiger]{gbresolve}} (if \code{taxonomy} is NULL)} 
}

\author{JM Eastman}

\examples{
sal=get(data(urodela))
print(head(sal$tax))

nphy=exemplar.phylo(sal$phy, strict.vals=NULL, sal$tax)
plot(nphy, type="fan", cex=0.15)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
