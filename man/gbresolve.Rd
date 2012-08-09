\name{gbresolve.phylo}
\alias{gbresolve.phylo}
\title{
NCBI taxonomy }
\description{
Working with NCBI taxonomy}
\usage{
	gbresolve.phylo(x, rank="phylum", within="", split=TRUE, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a phylogenetic tree of class 'phylo' or character string (taxonomic label)}
  \item{rank}{if non-NULL, a character vector specifying the columns of \code{taxonomy} to be used in relabeling}
  \item{within}{a character string representing a lineage within which to resolve the tips in 'x'}
  \item{split}{Boolean (default is \code{TRUE}); whether to use the first epithet encountered in tips labels of 'x'} 
  \item{...}{additional arguments to be passed to \code{load_taxdump}}
}

\author{JM Eastman}

\examples{
gbresolve("Ambystoma_laterale")
sal=get(data(urodela))
x=gbresolve(sal$phy)
plot(x$phy, show.node=TRUE, type="f", show.tip=FALSE, cex=0.3)
gbcontain("caudata", rank="genus")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
