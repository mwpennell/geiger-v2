\name{nh.test}
\alias{nh.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
using the Freckleton and Harvey node-height test
}
\description{
Fits a linear model between the absolute magnitude of the standardized independent contrasts and the height above the root of the node at which they were being compared to identify early bursts of trait evolution.
}
\usage{
nh.test(phy, d, regression.type, log = TRUE, rlm.maxit = 20 , show.plot = TRUE, \dots)
}
\arguments{

 \item{phy}{
A time calibrated phylogeny in "phylo" format
}
  \item{d}{
A named vector of trait values. 
}
  \item{regression.type}{
The type of regression to be used. Specify \code{regression.type="lm"} to fit a standard linear model or \code{regression.type="rlm"} to fit a robust regression model 
}
 
  \item{log}{
Whether the data should be logged or not.
}
  \item{rlm.maxit}{
The maximum number of iterations to fit the robust regression model. This is ignored if \code{regression.type="lm"}. 
}
  \item{show.plot}{
Binary argument indicating whether plot should be made.
}

  \item{...}{
Additional arguments passed to plot
}


}

\value{

Function returns a \code{\link[stats]{lm}} or \code{\link[MASS]{rlm}} object and outputs a plot if \code{show.plot=TRUE}
}
\references{
Slater GJ and MW Pennell (in press) Robust regression and posterior predictive simulation increase power to detect early bursts of trait evolution. \emph{Systematic Biology}.

Freckleton RP and PH Harvey (2006) Detecting non-brownian evolution in adaptive radiations. \emph{PLoS Biology} 4:e373.

}

\author{
Graham Slater
}

\seealso{
\code{\link[geiger]{pp.mcmc}}, \code{\link[ape]{pic}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\donttest{
data(whales)

tmp <- treedata(whales$phy, whales$dat[,1])

phy <- tmp$phy
dat <- tmp$data[,1]

nh.test(phy, dat, regression.type="lm", show.plot=TRUE)
}
}

