\name{mecca}
\alias{mecca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
running a MECCA analysis
}
\description{
Runs MECCA's hybrid ABC-MCMC algorithm to jointly estimate diversification rates and trait evolution from incompletely sampled comparative data. Many of the arguments taken by this function are the same as those in calibrateMecca().
}
\usage{
mecca(phy, richness, cladeMean, cladeVariance, model = c("BM", "Trend", "twoRate"),
   prior.list = list(priorSigma = c(-4.961845, 4.247066), priorMean = c(-10, 10)),
   start = start, Ngens = 10000, printFreq = 100, sigmaPriorType = "uniform",
   rootPriorType = "uniform", SigmaBounds = c(-4.961845, 4.247066), hotclade = NULL,
   divPropWidth = 0.1, scale = 1, divSampleFreq = 0, BoxCox = TRUE, outputName = "mecca")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
time-calibrated phylogenetic tree of class 'phylo'
}
  \item{richness}{
named vector of species richnesses corresponding to tips in the tree (if a tip does not represent a higher level taxon, its richness should be 1)
}
  \item{cladeMean}{
named vector of trait means (all tips in the tree must be represented)
}
  \item{cladeVariance}{
named vector of trait variances (all tips in the tree must be represented; if only one taxon is present, use 0 for the variance)
}
  \item{model}{
model of trait evolution to be used -- options currently implemented are "BM" = Brownian Motion, "Trend" = Brownian moton with a trend, and "twoRate" = two Brownian rate model (see \code{hotclade} below).
}
  \item{prior.list}{
a list containing prior distribution parameters (default values used if this argument is empty)
}
  \item{start}{
ouput of \code{\link[geiger]{startingpt.mecca}}
}
  \item{Ngens}{
number of generations to run MECCA
}
  \item{printFreq}{
frequency of printing the acceptance rate
}
  \item{sigmaPriorType}{
type of prior distribution on the Brownian rate parameter (currently either "uniform" or "normal")
}
  \item{rootPriorType}{
type of prior distribution on the root state rate parameter (currently "uniform" is available)
}
  \item{SigmaBounds}{
bounds for sigma (default values correspond to a wide range taken from Harmon et al. 2010)
}
  \item{hotclade}{
if a two-rate model is to be fit, this specifies which clade takes the second rate -- two names should be specified in a vector; either the two tip names that span the clade of interest, or the name of a terminal/internal edge and NULL if only one branch takes the second rate
}
  \item{divPropWidth}{
proposal width for the diversification MCMC (default value of 0.1 seems to work well)
}
  \item{scale}{
a numeric value by which the proposal width for trait evolution parameters will be multiplied (a value of 2 seems to work well, but this should be adjusted for each individual dataset)
}
  \item{divSampleFreq}{
whether new trees are simulated at every step -- the default (0) is yes; if a non-zero value is given, this will determine the frequency (every n steps) with which new tip trees are simulated
}
  \item{BoxCox}{
whether summaries are BOX-COX standardized -- default (1) is yes and is recommended; this should always be consistent with the calibration step
}
  \item{outputName}{
name stem for output file names
}
}
\details{
The output files produced are formatted to be used with the C++ Program ABCtoolbox (Wegmann et al. 2010), which produces adjusted posterior distributions and can perform model selection without likelihoods.
}
\value{
MECCA does not store any output in memory. Instead, five output files are generated to the current working directory. This files are fully compatible with ABCtoolbox (Wegmann et al. 2011). The first file (outputname_bdSimFile.tx) will output the posterior sample for diversification parameters. The second file (outputname_bmSimFile.txt) ouputs the sampled trait evolution parameters and their associated raw summary statistics while outputname_ObsFile.txt gives the observed summaries.For ABC toolbox though, it will often be more efficient to use pls-transformed versions of the observed and simulated summary statistics. These are available in outputname_distObs.txt, and outputname_distSimFile.txt.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Slater GJ, LJ Harmon, D Wegmann, P Joyce, LJ Revell, and ME Alfaro. 2012. Fitting models of continuous trait evolution to incompletely sampled comparative data using approximate Bayesian computation. \emph{Evolution} 66:752-762.
}

\author{
Graham Slater, Luke Harmon, Daniel Wegmann
}
\note{
The numbers printed to the screen during the run give the current generation, acceptance rate for diversification parameters, acceptance rate for trait evolutionary rate parameters and acceptance rate for root state parameters, respectively
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\dontrun{
data(carnivores)
phy <- carnivores$phy
data <- carnivores$dat
richness <- data[,1]
names(richness) <- rownames(data)

priors <- list(priorSigma = c(-4.5, 4.5), priorMean = c(-5, 2))

## CALIBRATION (far too short for a real analysis)
Cal <- calibrate.mecca(phy, richness, model = "BM", prior.list = priors, Ncalibrations = 1000)

params <- Cal$trait[, c(1,2)] ## extract the calibration BM parameters
stats <- Cal$trait[, -c(1,2)] ## extract the calibration summary stats

## now we run pls, determining combinations of summaries that explain variation in our parameters
## For BM, 2 components is sufficient. For more complex models, more componenets will be required.
require(pls)
myPlsr<-pls::plsr(as.matrix(params) ~ as.matrix(stats), scale=F, ncomp = 2)

plot(RMSEP(myPlsr)) ## Look at Root Mean Square error plots

summary(myPlsr) ## take a look at % variance explained

plsdat <- myPlsr$loadings

## extract means and variances from the carnivore data ##

cladeMean<-data[,2]
names(cladeMean)<-rownames(data)
cladeVariance<-data[,3]
names(cladeVariance)<-rownames(data)

## STARTING POINT
## And now we can compute starting values for the ABC-MCMC
start <- startingpt.mecca(Cal, phy, cladeMean, cladeVariance,
tolerance = 0.05, plsdat, BoxCox = TRUE)


## MECCA (far too short for a real analysis)
mecca(phy, richness, cladeMean, cladeVariance, model = "BM", prior.list = priors, start = start,
   Ngens = 1000, printFreq = 100, sigmaPriorType = "uniform", rootPriorType = "uniform",
   SigmaBounds = c(-4.5, 4.5), divPropWidth = 0.1, scale = 2, divSampleFreq = 0, BoxCox = TRUE,
   outputName ="MeccaBM.txt")

## PASTE UNCOMMENTED FOLLOWING LINE TO DROP FILES CREATED BY MECCA
 # unlink(dir(pattern=paste(r)),recursive=TRUE)

}
}
