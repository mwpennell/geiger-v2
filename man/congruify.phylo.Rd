\name{congruify.phylo}
\alias{congruify.phylo}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ultrametricization of trees from a supplied timetree}
\description{
automagically generating secondary calibrations}
\usage{
congruify.phylo(reference, target, taxonomy = NULL, tol = 0,
  scale=c(NA, "PATHd8", "treePL"), ncores=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reference}{an ultrametric tree used to time-scale the \code{target}}
  \item{target}{a phylogram that is sought to be ultrametricized based on the \code{reference} phylogeny}
  \item{taxonomy}{a linkage table between tips of the phylogeny and clades represented in the tree; rownames of 'taxonomy' should be tips found in the phylogeny}
  \item{tol}{branching time in \code{reference} above which secondary constraints will be applied to \code{target}}
  \item{scale}{\code{NA}, \code{"PATHd8"} or \code{"treePL"} (if \code{PATHd8} or \code{"treePL"} are available in the R \code{PATH})}
  \item{ncores}{number of cores to be used}
}
\details{
	This function uses the \code{reference} to inform secondary calibrations in the \code{target}. The primary output is a table of 'congruent' nodes between the \code{reference}
	and \code{target} with associated dates extracted from corresponding nodes in the \code{reference}.

	If multiple trees are supplied as the \code{reference}, a 'congruification' table is generated for each.

	If \code{scale="PATHd8"}, the \code{target} will be smoothed by \code{PATHd8} using the "d8 calculation" (see \url{http://www2.math.su.se/PATHd8/PATHd8manual.pdf}).
	This scaling method requires that \code{PATHd8} is available on the user's \code{PATH} variable that can be accessed by \code{\link[base]{Sys.getenv}("PATH")}.

  If \code{scale="treePL"}, the \code{target} will be smoothed by \code{treePL}.
  This scaling method requires that \code{treePL} is available on the user's \code{PATH} variable that can be accessed by \code{\link[base]{Sys.getenv}("PATH")}.
}

\references{
Eastman JM, LJ Harmon, and DC Tank. 2013. Congruification: support for time scaling large phylogenetic trees. \emph{Methods in Ecology and Evolution}, in press.
}

\author{JM Eastman}

\examples{
\dontrun{
sal=get(data(caudata))
res=congruify.phylo(sal$fam, sal$phy, sal$tax, tol=0, scale=NA, ncores=2)
print(res$calibrations)
plot(ladderize(sal$phy,right=FALSE), cex=0.35, type="fan", label.offset=2.5)
plot(ladderize(sal$fam,right=FALSE), cex=0.5, type="fan", label.offset=2.5, no.margin=FALSE)

# if you have PATHd8 installed you can also run
# res=congruify.phylo(sal$fam, sal$phy, sal$tax, tol=0, scale=\"PATHd8\")
# print(res)
}
}
\keyword{graphs}
\keyword{manip}
