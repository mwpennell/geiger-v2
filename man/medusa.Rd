\name{runMedusa}
\alias{runMedusa}
\alias{summaryMedusa}
\title{MEDUSA: Modeling evolutionary diversification using stepwise AIC}
\description{
	Fits a series of diversification models to phylogenetic and taxonomic data
	using stepwise AIC
}
\usage{
runMedusa(phy, richness, estimateExtinction=T, modelLimit=20, cutAtStem=T, startR = 0.05, startE = 0.5, ...)
summaryMedusa(phy, richness, out, cutoff=4, plotTree=T, useCorrection=F, cutAtStem=T) 
}
\arguments{
  \item{phy}{Phylogenetic tree in phylo format}
  \item{richness}{A table with taxon names in column 1 and extant diversities in column 2}
  \item{estimateExtinction}{If false, fits a series of pure-birth models; otherwise birth-death 
  		models are used}
  \item{modelLimit}{Maximum number of rate shifts to place on the tree}
  \item{cutAtStem}{When rate shifts are placed on branches in Medusa, they could be 
  		placed at either the beginning or the end of the stem branch. cutAtStem=T cuts at the 
  		beginning of the stem branch, while cutAtStem=F cuts at the end.}
  \item{...}{Additional arguments for nlm}
  \item{out}{Output from medusa, sent to summaryMedusa to compile results}
  \item{cutoff}{Cutoff value for differences in AIC scores when comparing models. More complex models
  				with an AIC score more than this number of units lower than simpler models are retained.}
  \item{plotTree}{Plot phylogenetic tree with branches under different models shown in different colors.}
  \item{useCorrection}{If true, use AICc instead of AIC}
  \item{startR}{Starting guess for net diversification rate. Try tweaking if Medusa crashes}
  \item{startE}{Starting guess for relative extinction rate. Try tweaking if Medusa crashes}
}
\details{
  MEDUSA fits a series of diversification models to a combination of phylogenetic and taxonomic data.
  The input is a phylogenetic tree with branch lengths proportional to time showing the relationship
  among clades, and the diversity for living species in all of those clades. All taxa missing from the 
  tree have to be assigned to one of the tip clades in the richness matrix.
  
  The algorithm first fits a Single diversification model to the entire dataset. Then, a series of breaks
  are added, so that different parts of the tree evolve with different parameter values (birth and or
  death rates). MEDUSA first compares all single-breakpoint models with the overall model, and selects the 
  best one. Then all possible two-breakpoint models are compared with the best single-breakpoint model, 
  and so on.
  
  summaryMedusa summarizes the results of MEDUSA. Here one must choose a cutoff value for improvement in AIC
  score necessary to accept more complex models.
	}
\value{
From medusa: a matrix where each row summarizes the next more complex model selected by the stepwise 
procedure. Columns record: the node where the next break was found, the likelihood, the number of 
parameters, the AIC score, and the AICc score (corrected for small sample size).

From summaryMedusa: a list summarizing the model selected using AIC (or AICc) scores and the cutoff
value. Each item in the list gives the likelihood for that part of the tree, two parameter estimates, 
r=b-d and epsilon=d/b, and the tip taxa included in that part of the tree.
	}
\references{ Alfaro, M. E., F. Santini, C. Brock, H. Alamillo, A. Dornburg, D. L. Rabosky, G. 
Carnevale, and L. J. Harmon. 2009. Nine exceptional radiations plus high turnover explain 
species diversity in jawed vertebrates. Proceedings of the National Academy of Sciences 
106:13410-13414.}
\author{Luke J. Harmon and Daniel Rabosky}
\seealso{

}
\examples{

data(geospiza)
attach(geospiza)

 

}
\keyword{arith}
