\name{phy.anova}
\alias{phy.anova}
\alias{phy.manova}
\title{phylogenetic ANOVA and MANOVA}
\description{
	computing phylogenetic ANOVA or MANOVA
}
\usage{
phy.anova(phy, data, group, nsim=1000)
phy.manova(phy, data, group, nsim=1000, test=c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"))

}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{data}{dependent variable(s) (should be continuous)}
  \item{group}{independent variable (should be a factor)}
  \item{nsim}{number of simulations for calculating p-value}
  \item{test}{the name of the test statistic to be used}
}
\details{

This function performs an ANOVA or MANOVA in a phylogenetic context.  First, the test statistic for ANOVA 
(one dependent variable) or MANOVA (more than one dependent variable) is calculated.  The null distribution 
of this test statistic is then obtained by simulating new sets of dependent variables on the phylogenetic
tree.  Simulations are run under a Brownian-motion model.  For ANOVA, the rate parameter is estimated from 
the average squared independent contrast; for MANOVA the simulations use an estimated variance-covariance (vcv) 
matrix from \code{\link[geiger]{ic.sigma}}.

For MANOVA, you can specify the test statistic for the summary table.  Wilks' statistic is most popular in the 
literature and for further details, see \code{\link[stats]{summary.manova}}.
}
\value{
	Standard ANOVA or MANOVA table and p-value based on simulations
	}
\references{
Garland Jr., T., A. W. Dickerman, C. M. Janis, and J. A. Jones. 1993. Phylogenetic analysis of covariance 
by computer simulation. \emph{Syst. Biol.} 42(3):265-292.
}
\author{ Luke J. Harmon }
\seealso{ \code{\link[stats]{anova}}; \code{\link[stats]{summary.manova}}; \code{\link[geiger]{vcv.phylo}} }
\examples{

geo=get(data(geospiza))
f<-as.factor(c(rep(0, 7), rep(1, 6)))

phy.manova(geo$phy, geo$dat, f)

x<-geo$dat[,"wingL"]

phy.anova(geo$phy, x, f) 

}

\keyword{ arith }
