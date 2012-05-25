\name{rc}
\alias{rc}
\title{ Relative cladogenesis test }
\description{
	Carries out relative cladogenesis test for all slices through the tree.	
}
\usage{
rc(phy, make.plot=TRUE, plot.bonf=FALSE, p.cutoff=0.05, cex=par("cex"))

}
\arguments{
  \item{phy}{ Phylogenetic tree }
  \item{make.plot}{ Plot tree with significant branches highlighted? }
  \item{plot.bonf}{ Use bonferroni correction when highlighting significant branches on the tree plot? (see below for more details) }
  \item{p.cutoff}{ Cutoff for significant p-values; only used when plotting }
  \item{cex}{ A numeric value giving the factor scaling of the tip and node labels (Character EXpansion). The default is to take 
  the current value from the graphical parameters. }
}
\details{
	A list of nodes is returned, along with the number of lineages alive just before that node, the 
	maximum number of descendents that any of those lineages has at the present day, a p-value for this
	observation under the null hypothesis of a birth-death process (that is, given the null, what is the 
	probability that one of these lineages had at least that many descendents), and the p-value after 
	Bonferroni correction (given that a total of n-1 comparisons are made).
	
	If a plot is made, asterisks will mark significantly diverse clades.  These asterisks appear just 
	to the right of the MRCA of the diverse clade.
	
	The Bonferroni correction used here is exceedingly conservative for a tree of any reasonable size. 
	I'm not sure that I recommend it, especially given the exploratory nature of this test and the non-
	independence of the comparisons. 
	
	One will often see significant results "trickle down" nodes in the tree - that is, if one clade is 
	expecially diverse, then one or more of its parent clades will also be diverse.  I think the most 
	parsimonious place to attribute this effect is to the most shallow significant branch - that is, 
	the branch closest to the tips (see Moore et al. 2004).

}
\value{
  Table of results with four columns: 
  Number of ancestors, Maximum descendents, p-value, Bonferroni-corrected p-value
}
\references{ Purvis, A., Nee, S. & Harvey, P. H. (1995) Proc. R. Soc. London Ser. B 260, pp. 329-333.
Moore, B.R., K.M.A. Chan, and M.J. Donoghue. 2004. Detecting diversification rate variation in supertrees. In O.R.P. Bininda-Emonds (ed.), Phylogenetic 
Supertrees: Combining Information to Reveal the Tree of Life, pp. 487-533. Kluwer Academic, Dordrecht, the Netherlands.
}
\author{ Luke J. Harmon }
\examples{

data(geospiza)
attach(geospiza)

rc(geospiza.tree)

}

\keyword{ arith }
