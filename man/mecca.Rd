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
mecca(phy, richness, cladeMean, cladeVariance, model = c("BM", "Trend", "twoRate"), prior.list = list(priorSigma = c(-4.961845, 4.247066), priorMean = c(-10, 10)), start = start, Ngens = 10000, printFreq = 100, sigmaPriorType = "uniform", rootPriorType = "uniform", SigmaBounds = c(-4.961845, 4.247066), hotclade = NULL, divPropWidth = 0.1, scale = 1, divSampleFreq = 0, BoxCox = TRUE, outputName = "mecca")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
A time calibrated phylogenetic tree of class "phylo"
}
  \item{richness}{
A named vector of species richnesses corresponding to tips in the tree; if a tip does not represent a higher level taxon, its richness should be 1. 
}
  \item{cladeMean}{
A names vector of trait means. All tips in the tree must be represented
}
  \item{cladeVariance}{
A names vector of trait variances. All tips in the tree must be represented. If only one taxon is present, use 0 for the variance
}
  \item{model}{
The model of trait evolution to be used. Options currently implemented are: "BM" = Brownian Motion, "Trend" = Brownian moton with a trend,"twoRate" = two Brownian rate model (see hotclade below).
}
  \item{prior.list}{
a list containing prior distribution parameters. If no values are specified, default values will be used.
}
  \item{start}{
The ouput of MeccaStartValues
}
  \item{Ngens}{
Number of generations to run MECCA for
}
  \item{printFreq}{
How frequently should Mecca print acceptance rates to the screen?
}
  \item{sigmaPriorType}{
The type of prior distribution on the Brownian rate parameter. Currently, uniform and normal are implemented.
}
  \item{rootPriorType}{
The type of prior distribution on the root state rate parameter. Currently, uniform is implemented.
}
  \item{SigmaBounds}{
If a normal prior is used for the rate, it can be bounded to reasonable values. Default values corresponding to a wide range taken from Harmon et al. (2010)
}
  \item{hotclade}{
If a two-rate model is to be fit, which clade takes the second rate. Two names should be specified in a vector; either the two tip names that span the clade of interest, or the name of a terminal/ internal edge and NULL if only one branch takes the second rate
}
  \item{divPropWidth}{
Proposal width for the diversification MCMC. The default value of 0.1 seems to work well.
}
  \item{scale}{
a numeric value by which the proposal width for trait evolution parameters will be multiplied. A value of 2 seems to work well, but this should be adjusted for each individual dataset
}
  \item{divSampleFreq}{
Should new trees be simulated at every step? The default (0) is yes. If a non-zero value is given, this will determine the frequency (every n steps) with which new tip trees are simulated. This option is designed for use with large datasets where simulating new trees will slow MECCA down considerably.
}
  \item{BoxCox}{
Should summaries be BOX-COX standardized. Default is yes and is recommended. This should always be consistent with the calibration step
}
  \item{outputName}{
Name stem for output file names
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
Slater, GJ, Harmon, LJ, Wegmann D, Joyce, P., Revell, LJ, Alfaro ME. in press. Evolution, Wegmann D et al. 2010. Bioinformatics 11:116.
}
\author{
Graham Slater, Luke Harmon, Daniel Wegmann
}
\note{
The numbers printed to the screen during the run give the current generation, acceptance rate for diversification parameters, acceptance rate for trait evolutionary rate parameters and acceptance rate for root state parameters, respectively
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\donttest{
data(carnivores);
phy <- carnivores$phy;
data <- carnivores$dat;

richness <- data[,1];
names(richness) <- rownames(data);

priors <- list(priorSigma = c(-4.5, 4.5), priorMean = c(-5, 2));

## CALIBRATION
Cal <- calibrate.mecca(phy, richness, model = "BM", prior.list = priors, Ncalibrations = 1000);

params <- Cal$trait[, c(1,2)]; ## extract the calibration BM parameters
stats <- Cal$trait[, -c(1,2)]; ## extract the calibration summary stats

## now we run pls, determining combinations of summaries that explain variation in our parameters

myPlsr<-plsr(as.matrix(params) ~ as.matrix(stats), scale=F, ncomp = 2); ## run PLS. For BM, 2 components is sufficient. For more complex models, more componenets will be required.

plot(RMSEP(myPlsr)); ## Look at Root Mean Square error plots

summary(myPlsr); ## take a look at % variance explained

plsdat <- myPlsr$loadings;

## extract means and variances from the carnivore data ##

cladeMean<-data[,2];
names(cladeMean)<-rownames(data);
cladeVariance<-data[,3];
names(cladeVariance)<-rownames(data);

## STARTING POINT
## And now we can compute starting values for the ABC-MCMC
start <- startingpt.mecca(Cal, phy, cladeMean, cladeVariance, tolerance = 0.05, plsdat, BoxCox = TRUE)


## MECCA
mecca(phy, richness, cladeMean, cladeVariance, model = "BM", prior.list = priors, start = start, Ngens = 10000, printFreq = 100, sigmaPriorType = "uniform", rootPriorType = "uniform", SigmaBounds = c(-4.5, 4.5), divPropWidth = 0.1, scale = 2, divSampleFreq = 0, BoxCox = TRUE, outputName ="MeccaBM.txt")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
