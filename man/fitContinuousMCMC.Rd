\name{fitContinuousMCMC}
\alias{fitContinuousMCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit models of continuous trait evolution to comparative data using MCMC
}
\description{
This function is essentially an MCMC version of the fitContinuous function that fits models using Maximum Likelihood. The main advantage of using MCMC is that informative prior distributions may be placed on node values when such information is available, for example from the fossil record. In such cases, model selection may be improved upon and differ from results using extant taxa only.
}
\usage{
fitContinuousMCMC(phy, d, model = c("BM", "Trend", "SSP", "ACDC.exp","ACDC.lin"),
Ngens = 1e+06, sampleFreq = 1000, printFreq = 1000,
propwidth = rep(0.1, 5), node.priors = NULL, root.prior = NULL,
acdc.prior = NULL, sample.node.states = T, outputName = "mcmc.output")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
A phylogenetic tree in phylo format
}
  \item{d}{
A named numeric vector of tip values
}
  \item{model}{
The type of model to be fitted. Options are: "BM" - Brownian motion, "Trend" = BM with a directional trend (this option can only be used with non ultrametric trees or when node priors are available for an ultrametric tree), "SSP" - the single stationary peak, or single-optimum OU model, "ACDC.exp" - Accelerating-Decelerating evolution where the change in rate is exponential with respect to time. The DC portion of this model corresponds to the EB model in Harmon et al. (2010), "ACDC.lin" - Accelerating-Decelerating evolution where the change in rate is linear with respect to time (see the supplementary information for Harmon et al. 2010 for details).
}
  \item{Ngens}{
The number of generations that the MCMC should be run for.
}
  \item{sampleFreq}{
The frequency with which likelihoods and parameters should be sampled from the chain. For large trees and/or if ancestral states are being sampled, this value should be large (at least equal to the number of internal nodes in the tree).
}
  \item{printFreq}{
The frequency with which progress should be printed to the screen.
}
  \item{propwidth}{
A numeric vector containing five values corresponding to proposal widths. The order is: 1 - root state proposal width, 2 - rate proposal width, 3 - scaler (alpha, rate change, trend parameter) proposal width, 4 - Theta (OU optimum) proposal width (not yet implemented), and 5- node state proposal width. 5 Values should be provided regardless of the model being fitted. Dummy or random values can be used for non-applicable parameters in such cases.
}
  \item{node.priors}{
A data frame of informative node priors. The default is node.priors=NULL, in which case uniformative priors will be assumed for all nodes. Alternatively, the user may provide a dataframe with five columns. The first two columns should specify the tips spanning the clade descended from the node on which an informative prior is to be placed. Columns 3 and 4 give parameters related to the prior while column 5 specifies the type of prior distribution. Current options, along with the associated parameters required in columns 3 and 4, are "normal" (mean and sd), "uniform" (min and max), and "exp" (offset and rate).
}
  \item{root.prior}{
A prior for the root value. The default is root.prior = NULL, in which case an uniformative root prior is assumed. Alternatively, the user may provide a vector containing 3 values: 1 & 2 are the parameters for the root prior parameter values, while the 3rd is the type of prior. Current prior types are the same as for node priors.
}
  \item{acdc.prior}{
For the ACDC models, it is  wise to bound the values of the rate change parameter as unreasonable proposed values can result in rounding issues and negatively infinite log-likelihoods. The default for this option is NULL. If an ACDC model is selected and no prior is given, a bounded, informative uniform prior will be assigned based on the root height of the tree. Alternatively, the user may provide bounds of their own. This may be particularly useful if the user is only interested in a subset of the parameter space area (eg., EB).
}
  \item{sample.node.states}{
Logical. If sample.node.states = TRUE, node states (ancestral state values) will be logged to a separate output file. If sample.node.states = FALSE, ancestral states will be treated as nuisance parameters and discarded after computation of likelihoods.
}
  \item{outputName}{
stem name for output file(s)
}
}
\details{
Unlike fitContinuous, which returns output as a list, fitContinuousMCMC outputs directly to text file(s). The model parameter file is designed to be compatatible with the popular Tracer software used in conjunction with BEAST. The output can easily be read back into R for further analysis. A particular advantage of writing directly to text file is that crashes etc do not result in the complete loss of long runs, and memory is saved.
}
\value{
If sample.node.states = FALSE, one text file is written to the current working directory. 
\item{outputName_model_params.txt}{a text file containing sampled prior, posterior, likelihood, and model parameter values from the MCMC run}
If sample.node.states = TRUE, a second text file is written
\item{outputName_nodestates.txt }{contains the posterior sample of ancestral states for internal nodes}
%% ...
}
\references{
Slater GJ, LJ Harmon, and ME Alfaro. 2012. Integrating fossils with molecular phylogenies improves inference of trait evolution. \emph{Evolution} 66:3931-3944.}
\author{
Graham Slater
}


\examples{
\dontrun{
data(caniformia)
phy <- caniformia$phy
d <- caniformia$dat
node.priors <- caniformia$node.priors
root.prior <- caniformia$root.prior

## as an example, we will run a very short (too short!) analysis,
##fitting BM and Trend to the caniform data

fitContinuousMCMC(phy, d, model = "BM", Ngens = 1000, sampleFreq =100,
printFreq = 100, node.priors = node.priors, root.prior = root.prior,
outputName ="BM_caniforms")

fitContinuousMCMC(phy, d, model = "Trend", Ngens = 1000, sampleFreq=100,
printFreq = 100, node.priors = node.priors, root.prior = root.prior,
outputName ="Trend_caniforms")

bm.res <- read.table("BM_caniforms_model_params.txt", header= TRUE)
head(bm.res)

trend.res <- read.table("Trend_caniforms_model_params.txt", header= TRUE)
head(trend.res)

### produce trace plots of logLk scores 

plot(bm.res$generation, bm.res$logLk, type = "l",
ylim = c(min(bm.res$logLk), max = max(trend.res$logLk)))

lines(trend.res$generation, trend.res$logLk, col = "red")
legend("bottomleft", c("bm", "trend"), lwd = 3, col = c("black", "red"))
}
}

