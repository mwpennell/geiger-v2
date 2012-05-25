\name{fitContinuous}
\alias{fitContinuous}
\title{ Model fitting for continuous data }
\description{
	Fits macroevolutionary models to phylogenetic trees
}
\usage{
fitContinuous(phy, data, data.names=NULL, model=c("BM", "OU", "lambda", "kappa", "delta", "EB", "white", "trend"), bounds=NULL,  meserr=NULL)
}
\arguments{
  \item{phy}{ object of type phylo }
  \item{data}{ Data vector (one trait) or matrix (multiple traits) }
  \item{data.names}{ Tip names for data vector that match tree species; ignored if data includes names}
  \item{model}{ Model to fit to comparative data; see below for options }
  \item{bounds}{ Range to constrain estimates; see below for details }
  \item{meserr}{ Measurement error for each trait for each tip species; can also be a single vector (if so, error is assumed to be equal for all traits) 
  					or a single number (error assumed equal for all traits across all species).}
}
\details{
This function fits various likelihood models for continuous character evolution.  The function 
returns parameter estimates, (approximate) confidence intervals based on the Hessian of the likelihood
function, and the likelihood.  Likelihood is maximized using the r function nlm.  This is a purely univariate function at this point - if multivariate data 
are sent to the function, it will carry out calculations for each character independently.  Possible models are
as follows:
 model="BM": Brownian motion 
 model="OU": Ornstein-Uhlenbeck; fits a model of random walk with a central tendency proportional to the parameter alpha.
 							Also called the "hansen" model in ouch, although the way the parameters are fit 
 							is slightly different here.
 model="lambda": Pagel's lambda; multiplies all internal branches of the tree by lambda, leaving tip branches as their original length.
 model="kappa": Pagel's kappa; raises all branch lengths to the power kappa.  As kappa approaches zero, the model becomes speciational.
 model="delta": Pagel's delta; raises all node depths to the power delta. If delta is less than one, evolution in concentrated early in the tree; delta > 1 concentrates evolution towards the tips.

 model="EB": Early burst, also called the ACDC model. Set by the eb parameter; fits a model where the rate of evolution increases or decreases exponentially through time, under the model
 					r(t) = ro * exp(r * t), where ro is the inital rate and r is the rate change parameter.  The actual parameter estimated, endRate, 
 					is the proportion of the initial rate represented by the end rate.
Bounds for the parameters in the likelihood search are sometimes necessary. This is because in some cases, the likelihood surface can have long flat ridges that 
cause the search to get "stuck."  This is particularly common, in my experience, for the OU model.  You can set bounds with the "bounds" parameter.  For example, 
to set bounds on alpha, use bounds=list(alpha=c(0.0001, 10)).  One can also set multiple bounds at once: bounds=list(alpha=c(0.0001, 10), beta=c(1, 100)). 
This function might work better than in previous versions of Geiger. 
model="white": White noise model; all species drawn from the same normal distribution (no phylogenetic signal).
model="trend": Brownian motion with a trend.

}
\value{
	Returns a matrix of parameter estimates, along with approximate standard errors and 95% confidence intervals. 
	Also returns the log-likelihood of the model (lnl).
}
\references{Lambda, kappa, and delta: Pagel, M. 1999. Inferring the historical patterns of biological evolution. Nature 401:877-884.
			OU: Butler, M.A. and A.A. King, 2004. Phylogenetic comparative analysis: a modeling approach for adaptive evolution. American Naturalist 164:683-695.
			Early burst: Paper in revision., L. J. Harmon et al. }
\author{Luke J. Harmon and Wendell Challenger}
\examples{

data(geospiza)
attach(geospiza)


#---- PRINT RESULTS
fitContinuous(geospiza.tree, geospiza.data)
 
#---- STORE RESULTS 
brownFit <-  fitContinuous(geospiza.tree, geospiza.data)

aic.brown<-numeric(5)
for(i in 1:5) aic.brown[i]<-brownFit[[i]]$aic

#----------------------------------------------------
#   PHYLOGENETIC SIGNAL: FIT LAMBDA
#----------------------------------------------------

lambdaFit<-fitContinuous(geospiza.tree, geospiza.data, model="lambda")

# Compare likelihoods:

d.lambda<-numeric(5)
for(i in 1:5) d.lambda[i]=2*(lambdaFit[[i]]$lnl-brownFit[[i]]$lnl)

# Calculate p values assuming chi-squared distribution with 1 d.f.
p.lambda=pchisq(d.lambda, 1, lower.tail=FALSE)

aic.lambda<-numeric(5)
for(i in 1:5) aic.lambda[i]<-lambdaFit[[i]]$aic

#----------------------------------------------------
#    TIME PROPORTIONALITY: DELTA 
#---------------------------------------------------

deltaFit<-fitContinuous(geospiza.tree, geospiza.data, model="delta")

# Compare likelihoods:

d.delta<-numeric(5)
for(i in 1:5) d.delta[i]=2*(deltaFit[[i]]$lnl-brownFit[[i]]$lnl)

# Calculate p values assuming chi-squared distribution with 1 d.f.
p.delta=pchisq(d.delta, 1, lower.tail=FALSE)

aic.delta<-numeric(5)
for(i in 1:5) aic.delta[i]<-deltaFit[[i]]$aic

#----------------------------------------------------
#   SPECIATIONAL MODEL: KAPPA 
#---------------------------------------------------

kappaFit<-fitContinuous(geospiza.tree, geospiza.data, model="kappa")

# Compare likelihoods:

d.kappa<-numeric(5)
for(i in 1:5) d.kappa[i]=2*(kappaFit[[i]]$lnl-brownFit[[i]]$lnl)

# Calculate p values assuming chi-squared distribution with 1 d.f.
p.kappa=pchisq(d.kappa, 1, lower.tail=FALSE)

aic.kappa<-numeric(5)
for(i in 1:5) aic.kappa[i]<-kappaFit[[i]]$aic

#----------------------------------------------------
#   OU MODEL: ALPHA 
#---------------------------------------------------

ouFit<-fitContinuous(geospiza.tree, geospiza.data, model="OU")

# Compare likelihoods:

d.ou<-numeric(5)
for(i in 1:5) d.ou[i]=2*(ouFit[[i]]$lnl-brownFit[[i]]$lnl)

# Calculate p values assuming chi-squared distribution with 1 d.f.
p.ou=pchisq(d.ou, 1, lower.tail=FALSE)

aic.ou<-numeric(5)
for(i in 1:5) aic.ou[i]<-ouFit[[i]]$aic

#----------------------------------------------------
#   EARLY BURST MODEL: R 
#---------------------------------------------------

ebFit<-fitContinuous(geospiza.tree, geospiza.data, model="EB")

# Compare likelihoods:

d.eb<-numeric(5)
for(i in 1:5) d.eb[i]=2*(ebFit[[i]]$lnl-brownFit[[i]]$lnl)

# Calculate p values assuming chi-squared distribution with 1 d.f.
p.eb=pchisq(d.eb, 1, lower.tail=FALSE)

aic.eb<-numeric(5)
for(i in 1:5) aic.eb[i]<-ebFit[[i]]$aic

#----------------------------------------------------
#   COMPARE ALL MODELS
#---------------------------------------------------

# One way: use likelihood ratio test to compare all models to Brownian model

d.all<-cbind(d.lambda, d.delta, d.kappa, d.ou, d.eb)
p.all<-cbind(p.lambda, p.delta, p.kappa, p.ou, p.eb)

cat("Trait\tlambda\tdelta\tkappa\tou\teb\n")

for(i in 1:5) {
	cat("Tr", i, "\t");
	for(j in 1:5) {
		cat(round(d.all[i,j],2));
		if(p.all[i,j]<0.05) cat("*");
		if(p.all[i,j]<0.01) cat("*");
		if(p.all[i,j]<0.001) cat("*");
		cat("\t");
	}
	cat("\n");
}

# Another way: use AIC

aic.all<-cbind(aic.brown, aic.lambda, aic.delta, aic.kappa, aic.ou, aic.eb)
foo<-function(x) x-x[which(x==min(x))]
daic<-t(apply(aic.all, 1, foo))

rownames(daic)<-colnames(geospiza.data)
colnames(daic)<-c("Brownian", "Lambda", "Delta", "Kappa", "OU", "EB")

cat("Table of delta-aic values; zero - best model\n")
print(daic, digits=2)



}

\keyword{ arith }
