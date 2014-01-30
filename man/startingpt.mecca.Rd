\name{startingpt.mecca}
\alias{startingpt.mecca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
starting values for MECCA
}
\description{
This function takes the output of calibrateMecca along with observed data and partial least squares loadings and outputs starting values and tuning parameters for the ABC-MCMC (Wegmann et al. 2009).
}
\usage{
startingpt.mecca(calibrationOutput, phy, cladeMean, cladeVariance,
                 tolerance = 0.01, plsComponents, BoxCox = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{calibrationOutput}{
The output from calibrateMecca
}
  \item{phy}{
A time calibrated phylogeny in "phylo" format
}
  \item{cladeMean}{
A named vector of trait means. All tips in the tree must be represented
}
  \item{cladeVariance}{
A names vector of trait variances. All tips in the tree must be represented. If only one taxon is present, use 0 for the variance
}
  \item{tolerance}{
The proportion of calibrations simulations that fall closest to the oberved data that will be retained to compute MECCA tuning parameters
}
  \item{plsComponents}{
a matrix of Partial Least Squares component loadings
}
  \item{BoxCox}{
Logical - Should Summary Statistics be Box-Cox transformed? Default is yes and is recommended
}
}
\details{
You will need to compute PLS loadings using the package "pls" prior to running this function. MECCA performs extremely poorly if summaries are not PLS transformed. If Bayes Factors are to be computed to perform model selection (Leuenberger and Wegmann 2010), raw summaries will need to be used in the post-sampling adjustment step. However, PLS transformed summaries can still be used in the acceptance/rejection step of the MCMC and can also be used to determine which simulations to retain
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{$tuning }{a matrix containing starting values and proposal widths for trait evolution parameters}
 \item{$startingBirth }{Starting value for speciation}
 \item{$startingDeath }{Starting value for extinction}
 \item{$dcrit }{the criticial distance - simulated data must fall within this distance of the observed data in order to be accepted}
 \item{$obsTraits }{the observed data}
 \item{$plsObserved }{Partial Least Squared transformed observed data}
 \item{$plsLoadings }{Partial Least Squared componenent loadings - these are used to transform data during the MCMC}
}
\references{
Slater GJ, Harmon LJ, Wegmann D, Joyce P, Revell LJ, Alfaro ME. in press Evolution, Leuenberger C, and Wegmann D. 2010. Genetics 184: 243-252., Wegmann D, Leuenberger C, Excoffier L. 2009. Genetics 182: 1207-1218.
}
\author{
Graham Slater, Luke Harmon, Daniel Wegman
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\donttest{
example(mecca)
}
}
