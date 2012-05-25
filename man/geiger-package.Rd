\name{geiger-package}
\alias{geiger-package}
\alias{geiger}
\docType{package}
\title{
GEIGER
}
\description{
 A package for macroevolutionary simulation and estimating parameters 
 related to diversification from comparative phylogenetic data.
}
\details{
\tabular{ll}{
Package: \tab geiger\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-01-19\cr
License: \tab GPL version 2 or greater?\cr
}

}
\author{
Luke J. Harmon, Jason Weir, Chad Brock, and Wendell Challenger

Maintainer: Luke Harmon <lukeh@uidaho.edu>
}
\examples{

# Disparity-through-time, as in Harmon et al. 2003
data(geospiza)
attach(geospiza)

dtt.full(geospiza.tree, geospiza.data)


# Simulation tester - verifies that simulations produce expected vcv structure
sims<-sim.char(geospiza.tree, as.matrix(1.0), n=100000)
m<-var(t(sims[,1,]))
m2<-vcv.phylo(geospiza.tree)
round(m-m2)
# Equal to zero because the vcv matrices are as expected.

#Test models of evolution
fitContinuous(geospiza.tree, geospiza.data)


# Reconstruct evolutionary vcv matrix
v<-cbind(c(4.0, -1.5),c(-1.5, 9.0))
n<-1000
sims<-sim.char(geospiza.tree, v, n)
r<-array(dim=c(2, 2, n))
for(i in 1:n)
	r[,,i]<-ic.sigma(geospiza.tree, sims[,,i])
v.recon<-apply(r, c(1,2), mean)
v-v.recon

# Discrete character evolution
q<-list(rbind(c(-.01, .01), c(.01, -.01)))
sims<-sim.char(geospiza.tree, q, model="discrete", n=100)

}
\keyword{arith}
\keyword{graphs}