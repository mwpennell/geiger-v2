\name{congruify.phylo}
\alias{congruify.phylo}
\alias{nodefind.phylo}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ultrametricization of trees from a supplied timetree}
\description{
Generates a table of secondary calibrations}
\usage{
congruify.phylo(stock, scion, taxonomy = NULL, tol = 0, method="PATHd8")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stock}{an ultrametric tree used to time-scale the \code{scion}}
  \item{scion}{a phylogram that is sought to be ultrametricized based on the \code{stock} phylogeny}
  \item{taxonomy}{a linkage table between tips of the phylogeny and clades represented in the tree; rownames of 'taxonomy' should be tips found in the phylogeny}
  \item{tol}{branching time in \code{stock} above which secondary constraints will be applied to \code{scion}}
  \item{method}{either \code{NA} or \code{"PATHd8"} (if \code{PATHd8} is available in the system \code{$PATH})}
}
\details{
	This function uses the \code{stock} to inform secondary calibrations in the \code{scion}.  The primary output is a table of 'congruent' nodes between the \code{stock} 
	and \code{scion} with associated dates given by the \code{stock}.  
	
	If multiple trees are supplied as the \code{stock}, a 'congruification' table is generated for each.  
	
	If \code{method="PATHd8"}, the \code{scion} will be smoothed by \code{PATHd8} using the "d8 calculation" (see \url{http://www2.math.su.se/PATHd8/PATHd8manual.pdf}).
}

\author{JM Eastman}

\examples{
sal=get(data(urodela))
res=congruify.phylo(sal$fam, sal$phy, sal$tax, tol=0, method=NA)
print(res$calibrations)
plot.phylo(ladderize(sal$phy,right=FALSE), cex=0.35, type="fan", label.offset=2.5)
plot.phylo(ladderize(sal$fam,right=FALSE), cex=0.5, type="fan", label.offset=2.5, no.margin=FALSE)

cat("\n\n\n*** If 'PATHd8' is installed ***, execute the following:\n\tres=congruify.phylo(sal$fam, sal$phy, sal$tax, tol=0, method=\"PATHd8\")\n\tprint(res)\n\tplot(res$phy, cex=0.35, type=\"fan\", label.offset=2.5)\n\n\tprint(max(branching.times(res$phy)))\n\tprint(max(branching.times(sal$fam)))\n\tprint(max(branching.times(sal$phy)))\n\n\n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
