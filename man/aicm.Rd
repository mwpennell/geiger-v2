\name{aicm}
\alias{aicm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Akaike's Information Criterion for MCMC samples (AICM)
}
\description{
Computes Akaike's Information Criterion for MCMC samples (AICM: see Raftery et al. 2007). Can be used to perform model selection using output from fitContinuousMCMC.
}
\usage{
aicm(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vector containing sampled likelihoods from the MCMC chain. Assumes that burn-in has been removed prior to computation of AICM score.
}
}
\details{
AICM is one way of comparing model fit using posterior likelihood samples. It has advantages over approaches such as thermodynamic integration in that it uses the chain output directly and thus has little added time cost. Furthermore, it has been shown to perform better than the harmonic mean estimator of the marginal likelihood. However, it is also less robust than stepping-stone or thermodynamic integration approaches and should be used with care.
}
\value{
%%  ~Describe the value returned
AICM - Akaike's Information Criterion for the posterior sample 
}
\references{
Raftery et al. 2007. Estimating the integrated likelihood via posterior simulation using the harmonic mean identity. In. Bernardo et al. (eds) Bayesian Statistics. Oxford University Press.
}
\author{
Graham Slater
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

## generate a random set of values from a normal distribution, 
## similar to a set of likelihood samples obtained via MCMC.

x <- rnorm(1000, -275, 2);
aicm(x);

}
