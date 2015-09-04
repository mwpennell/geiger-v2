\name{gbresolve}
\alias{gbresolve.phylo}
\alias{gbresolve}
\alias{gbcontain}

\title{ NCBI taxonomy }
\description{
working with NCBI taxonomy}
\usage{
	gbresolve(x, rank="phylum", within="", ncores=1, ...)
	gbcontain(x, rank="species", within="", ncores=1,...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a phylogenetic tree of class 'phylo' (\code{gbresolve.phylo}), or a string vector (\code{gbresolve}), or a single string (\code{gbcontain})}
  \item{rank}{a Linnaean rank to restrict taxonomic resolution}
  \item{within}{a character string representing a group within which to resolve query}
  \item{ncores}{number of cores to use.}
	\item{...}{additional arguments to be passed to the taxdump constructor (\code{\link[ncbit]{ncbit}})}
}

\details{
	The functions access the NCBI taxonomy resource (\url{http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html}; see also \code{\link[ncbit]{ncbit}}).
	\code{gbresolve} resolves the taxonomic hierarchy for queried taxa up to the given \code{rank} (or between the given \code{rank}s if two are given),
    and \code{gbcontain} resolves all taxa found within a given queried group and occurring at a specified \code{rank}.  The \code{rank} must be found
    within the object \code{Linnaean} (see \bold{Examples}). The argument \code{within} can restrict the group within which to conduct the search
    (see \bold{Examples}).

	The local copy of the taxonomy resource (accessible in with \code{data(ncbi)} from \pkg{ncbit}) can be updated with a call to \code{ncbit(update=TRUE)}.

	Setting the \code{ncores} argument to \code{NULL} will use all available cores.
}

\author{JM Eastman}

\examples{
\donttest{
## call up NCBI taxonomy
ncbi=ncbit::ncbit(update=FALSE)

## possible ranks
print(Linnaean)

## resolve taxa
gbresolve(c("Ambystoma_laterale", "Dicamptodon_copei"))
gbresolve("Andrias japonicus")


## resolve taxa found in tree
sal=get(data(caudata))
x=gbresolve(sal$phy, rank=c("genus", "order"))
plot(x$phy, show.node=TRUE, type="f", show.tip=FALSE, cex=0.3)

## find all genera within salamanders
gbcontain("caudata", rank="genus")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data}
\keyword{graphs}
