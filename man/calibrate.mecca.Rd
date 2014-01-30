\name{calibrate.mecca}
\alias{calibrate.mecca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calibrating MECCA
}
\description{
Runs a user defined number of calibration simulations (Wegmann et al. 2009) to generate tuning parameters MECCA's ABC-MCMC. The function takes a large number of arguments but in most cases these can remain at default settings. 
}
\usage{
calibrate.mecca(phy, richness, model = c("BM", "Trend", "twoRate"), 
	prior.list = list(priorSigma = c(-4.961845, 4.247066), priorMean = c(-10, 10)), 
	Ncalibrations = 10000, sigmaPriorType = "uniform", rootPriorType = "uniform", 
	divSampleFreq = 0, initPars = "ML", propWidth = 0.1,
        SigmaBounds = c(-4.961845, 4.247066), 
	hotclade = NULL, BOXCOX = TRUE, optimRange = c(-1000, 10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
A time calibrated phylogenetic tree of class "phylo"
}
  \item{richness}{
A named vector of species richnesses corresponding to tips in the tree; if a tip does not represent a higher level taxon, its richness should be 1. 
}
  \item{model}{
The model of trait evolution to be used. Options currently implemented are: "BM" = Brownian Motion, "Trend" = Brownian moton with a trend, "twoRate" = two Brownian rate model (see hotclade below).
}
  \item{prior.list}{
A namedlist containing prior distribution parameters. If no values are specified, default values will be used. Default values for the BM rate are taken from the range reported in Harmon et al. 2010.
}
  \item{Ncalibrations}{
The number of calibration simulations to use. No fewer than 10,000 are recommended.
}
  \item{sigmaPriorType}{
The type of prior distribution on the Brownian rate parameter. Currently, uniform and normal are implemented.
}
  \item{rootPriorType}{
The type of prior distribution on the root state rate parameter. Currently, uniform is implemented.
}
  \item{divSampleFreq}{
Should new trees be simulated at every step? The default (0) is yes. If a non-zero value is given, this will determine the frequency (every n steps) with which new tip trees are simulated. This option is designed for use with large datasets where simulating new trees will slow MECCA down considerably.
}
  \item{initPars}{
Should Maximum Likelihood ("ML") values be used to start the diversification rate sampler? if not, two numeric values should be specified for speciation and extinction rates.
}
  \item{propWidth}{
Proposal width for the diversification MCMC. The default value of 0.1 seems to work well.
}
  \item{SigmaBounds}{
If a normal prior is used for the trait evolution rate, it can be bounded with reasonable values to ensure that the chain does not get stuck in areas of low likelihood. Default values correspond to the range of rates reported in Harmon et al. (2010).
}
  \item{hotclade}{
If a two-rate model is to be fit, which clade takes the second rate? There are two ways to do this - either a clade is designated or a single edge. For a clade, two tip names that span the clade of interest should be specified in a vector. If a single edge is to be used pass a two-item vector, where the first item is a single terminal name (for a tip branch) or the number of a node (for an internal edge) and the second item is NULL
}
  \item{BOXCOX}{
Will summaries be BOX-COX standardized? Default is yes and is recommended. Raw summaries will be returned, along with parameters for Box-Cox transformation.
}
  \item{optimRange}{
Optimization range for the lambda parameter for Box-Cox. Default params should work in most cases.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{$diversification }{A 3-column matrix containing 1) birth-rates and 2) death-rates sampled by the MCMC, plus 3) their associated likelihoods}
 \item{$trait }{A matrix of trait evolution parameters sampled from their priors and summary statistics from the simulations. If BM is used, columns 1 and 2 are the rate and root state parameter samples. If Trend or twoRate are the models, the first two columns are the same but the third contains the additional rate or trend parameter}
 \item{$stdz, $lambda, $GM, $BoxCox }{If BOXCOX = T, these items will also be output to BoxCox transform simulations in the MCMC}

}
\references{
Slater, GJ, Harmon, LJ, Wegmann D, Joyce, P., Revell, LJ, Alfaro ME. in press. Evolution, also Wegmann D, Leuenberger C, Excoffier L. 2009. Genetics 182:1207-1218.
Harmon LJ et al. 2010. Early Bursts of body size and shape evolution are rare in comparative data. Evolution 64: 2385 - 2396}
\author{
Graham J Slater, Luke Harmon, Dan Wegmann
}
\note{
Tree simulation is done using modified code from Tanja Stadler's TreeSim package. Trait evolution is done using modified code from Liam Revell's phytools package.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\donttest{
example(mecca)
}
}

