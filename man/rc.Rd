\name{rc}
\alias{rc}
\title{ relative cladogenesis test }
\description{
	conducting the relative cladogenesis test for all slices through a tree	
}
\usage{
rc(phy, plot=TRUE, bonf=FALSE, p.cutoff=0.05, cex=par("cex"))

}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{plot}{ whether to plot tree with significant branches highlighted}
  \item{bonf}{ whether to use bonferroni correction when highlighting significant branches on the tree (see \bold{Details})}
  \item{p.cutoff}{cutoff for significant p-values (used only when plotting)}
  \item{cex}{ a numeric value giving the factor scaling of the tip and node labels (see \code{\link[graphics]{par}})}
}
\details{
	A list of nodes is returned, along with the number of lineages alive just before that node, the 
	maximum number of descendents that any of those lineages has at the present day, a p-value for this
	observation under the null hypothesis of a birth-death process (that is, given the null, what is the 
	probability that one of these lineages had at least that many descendents), and the p-value after 
	Bonferroni correction (given that a total of n-1 comparisons are made).
	
	If a plot is made, asterisks will mark significantly diverse clades.  These asterisks appear just 
	to the right of the MRCA of the diverse clade.
	
	The Bonferroni correction used here is exceedingly conservative for a tree of any reasonable size (and not 
	necessarily recommended, especially given the exploratory nature of this test and the non-
	independence of the comparisons). 
	
	One will often see significant results "trickle down" nodes in the tree - that is, if one clade is 
	expecially diverse, then one or more of its parent clades will also be diverse.  The most 
	parsimonious place to attribute this effect is to the most shallow significant branch - that is, 
	the branch closest to the tips (see Moore et al. 2004).

}
\value{
  Table of results with four columns: 
  Number of ancestors, Maximum descendents, p-value, Bonferroni-corrected p-value
}
\references{ 
Purvis, A., Nee, S. & Harvey, P. H. 1995. \emph{Proc. R. Soc. London Ser. B} 260:329-333.

Moore, B.R., K.M.A. Chan, and M.J. Donoghue. 2004. Detecting diversification rate variation in supertrees. In O.R.P. Bininda-Emonds (ed.), \emph{Phylogenetic 
Supertrees: Combining Information to Reveal the Tree of Life}, pp. 487-533. Kluwer Academic, Netherlands:Dordrecht.
}
\author{ Luke J. Harmon }
\examples{

geo=get(data(geospiza))

rc(geo$phy)

}

\keyword{ arith }
