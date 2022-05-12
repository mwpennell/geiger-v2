\name{aov.phylo}
\alias{aov.phylo}
\title{phylogenetic ANOVA and MANOVA}
\description{
	computing phylogenetic ANOVA or MANOVA
}
\usage{
aov.phylo(formula, phy, nsim = 1000, 
    test = c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula specifying the model (see \bold{Examples})}
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{nsim}{number of simulations to run}
  \item{test}{test statistic to apply if MANOVA}
  \item{\dots}{additional arguments to be passed to \code{print.anova}}
}
\details{
This function performs an ANOVA or MANOVA in a phylogenetic context.  First, the test statistic for ANOVA 
(one dependent variable) or MANOVA (more than one dependent variable) is calculated.  The null distribution 
of this test statistic is then obtained by simulating new sets of dependent variables on the phylogenetic
tree.  Simulations are run under a Brownian-motion model.  For ANOVA, the rate parameter is estimated from 
the average squared independent contrast; for MANOVA the simulations use an estimated variance-covariance
matrix from \code{\link[geiger]{ratematrix}}.

For MANOVA, you can specify the test statistic for the summary table.  Wilks' statistic is most popular in the 
literature and for further details, see \code{\link[stats]{summary.manova}}.
}
\value{
The function prints (and returns) a standard ANOVA or MANOVA table and p-value based on simulations (from the \code{Pr(phy)} column). For convenience, the summary table is included as an attribute of the returned object (see \bold{Examples}).
}
\references{
Garland T Jr, AW Dickerman, CM Janis, and JA Jones. 1993. Phylogenetic analysis of covariance by computer simulation. \emph{Systematic Biology} 42(3):265-292.
}
\author{ JM Eastman and LJ Harmon }
\seealso{ \code{\link[stats]{aov}}; \code{\link[stats]{anova}}; \code{\link[stats]{summary.manova}}; \code{\link[geiger]{ratematrix}} }
\examples{
\dontrun{
geo=get(data(geospiza))
dat=geo$dat
d1=dat[,1]
grp<-as.factor(c(rep(0, 7), rep(1, 6)))
names(grp)=rownames(dat)

## MANOVA
x=aov.phylo(dat~grp, geo$phy, nsim=50, test="Wilks")
print(attributes(x)$summary) # summary table

## ANOVA
x1=aov.phylo(d1~grp, geo$phy, nsim=50)
}
}

\keyword{univar}
\keyword{multivariate}


