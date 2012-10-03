\name{gbresolve}
\alias{gbresolve.phylo}
\alias{gbresolve}
\alias{gbcontain}
\alias{Linnaean}

\title{ NCBI taxonomy }
\description{
working with NCBI taxonomy}
\usage{
	gbresolve(x, rank="phylum", within="", split=FALSE, ...)
	gbcontain(x, rank="species", within="", ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a phylogenetic tree of class 'phylo' or a string vector (\code{gbresolve}) or a single string (\code{gbcontain})}
  \item{rank}{a Linnaean rank to restrict taxonomic resolution}
  \item{within}{a character string representing a group within which to resolve query}
  \item{split}{whether to use the first epithet encountered in tips labels of 'x'} 
  \item{...}{additional arguments to be passed to the taxdump constructor (\code{geiger:::.ncbi})}
}

\details{
	The functions access the NCBI taxonomy resource (\url{http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html}; see also \code{\link[geiger]{ncbi}}). 
	\code{gbresolve} resolves the taxonomic hierarchy for queried taxa up to the given \code{rank}, and \code{gbcontain} resolves all taxa found within 
	a given queried group and occurring at a specified \code{rank}.  The \code{rank} must be found within the object \code{Linnaean} (see \bold{Examples}). 
	The argument \code{split=TRUE} starts resolution from the first epithet in the queried names (e.g., the genus epithet \emph{Ambystoma} in \emph{Ambystoma laterale}). 
	The argument \code{within} can restrict the group within which to conduct the search (see \bold{Examples}).  
	
	The local copy of the taxonomy resource (accessible in with \code{data(ncbi)}) can be updated with a call to \code{geiger:::.ncbi(update=TRUE)}.  
}

\author{JM Eastman}

\examples{
## possible ranks
print(Linnaean)

## resolve taxa
gbresolve(c("Ambystoma_laterale", "Dicamptodon_copei"))
gbresolve("Andrias japonicus")


## resolve taxa found in tree
sal=get(data(caudata))
x=gbresolve(sal$phy)
plot(x$phy, show.node=TRUE, type="f", show.tip=FALSE, cex=0.3)

## find all genera within salamanders
gbcontain("caudata", rank="genus")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line