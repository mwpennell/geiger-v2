\name{gbresolve}
\alias{gbresolve}
\title{
NCBI taxonomy }
\description{
Working with NCBI taxonomy}
\usage{
	gbresolve(x, rank="phylum", within="", update=FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a phylogenetic tree of class 'phylo' or character string (taxonomic label)}
  \item{rank}{if non-NULL, a character vector specifying the columns of \code{taxonomy} to be used in relabeling}
  \item{within}{a ('matrix') linkage table between tips of the phylogeny and clades represented in the tree; rownames of 'taxonomy' should be tips found in the phylogeny}
  \item{update}{Boolean (default is \code{FALSE}); whether to rebuild a local copy of the NCBI taxonomy database} 
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
