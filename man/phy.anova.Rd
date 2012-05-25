\name{phy.anova}
\alias{phy.anova}
\alias{phy.manova}
\title{Phylogenetic ANOVA and MANOVA}
\description{
	Calculates ANOVA or MANOVA, and returns p-value based on Brownian motion simulations
}
\usage{
phy.anova(phy, data, group, data.names=NULL, nsim=1000)
phy.manova(phy, data, group, data.names=NULL, nsim=1000, test="Wilks")

}
\arguments{
  \item{phy}{ an object of class 'phylo'}
  \item{data}{ Dependent variable(s) - should be continuous }
  \item{group}{ Independent variable - should be a factor }
  \item{data.names}{ Taxon names in order, corresponding to x and group; not needed if data has names stored as 'names' or 'rownames'}
  \item{nsim}{ Number of simulations for calculating p-value; default = 1000 }
  \item{test}{ If you're using a MANOVA, the name of the test statistic to be used. 
  				Partial matching is used so the name can be abbreviated. 
  				Options are "Pillai", "Wilks", "Hotelling-Lawley", and "Roy"}
}
\details{

This function allows an ANOVA or MANOVA in a phylogenetic context.  First, the test statistic for ANOVA 
(one dependent variable) or MANOVA (more than one dependent variable) is calculated.  The null distribution 
of this test statistic is then obtained by simulating new sets of dependent variables on the phylogenetic
tree.  Simulations are run under a Brownian motion model.  For ANOVA, the rate parameter is estimated from 
the average squared independent contrast; for MANOVA the simulations use an estimated variance-covariance (vcv) 
matrix from the GEIGER function ic.sigma.

For MANOVA, you can specify the test statistic for the summary table.  Wilks' statistic is most popular in the 
literature; for more details see the summary.manova help page.
}
\value{
	Standard ANOVA or MANOVA table and p-value based on simulations
	}
\references{Garland Jr., T., A. W. Dickerman, C. M. Janis, and J. A. Jones. 1993. Phylogenetic analysis of covariance 
by computer simulation. Syst. Biol. 42(3):265-292.}
\author{ Luke J. Harmon }
\seealso{ anova, summary.manova, ic.sigma }
\examples{

data(geospiza)
attach(geospiza)


f<-as.factor(c(rep(0, 7), rep(1, 6)))


phy.manova(geospiza.tree, geospiza.data, f)

x<-geospiza.data[,1]

phy.anova(geospiza.tree, x, f, data.names=rownames(geospiza.data)) 

}

\keyword{ arith }
