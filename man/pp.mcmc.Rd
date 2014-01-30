\name{pp.mcmc}
\alias{pp.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
using posterior predictive MCMC for modeling quantitative trait evolution
}
\description{
performs posterior predictive checks for models of quantitative trait evolution. At present, only BM, EB, and clade.shift models are implemented
}
\usage{
pp.mcmc(phy, d, Ngens = 1000000, sampleFreq = 1000, printFreq = 1000,
        prop.width = 1, model = "BM", eb.type = "exponential",
        clade = NULL, rlm.maxit = 20)
}
\arguments{

 \item{phy}{
A time calibrated phylogeny in "phylo" format
}
  \item{d}{
A named vector or dataframe of trait values. 
}
  \item{Ngens}{
Number of generations that the posterior predictive MCMC will run for. Default is 1 million generations 
}
 
  \item{sampleFreq}{
The frequency with which model parameters will be sampled from the chain and simulations run. Default is every 1000 generations
}
  \item{printFreq}{
The frequency with which the current number of generations and acceptance rates will be printed to screen. Default is every 1000 generations
}
  \item{prop.width}{
The width of the sliding window proposal distribution for ln(Sigmasq) and, if applicable, the exponential change parameter for EB. The width for the EB parameter is obtained by dividing by 10. Default proposal width is 1.
}

  \item{model}{
The model to fit and simulate under. Default is Brownian motion (BM). Other options are early burst (EB) or an edge shift model (edge.shift) where the rate is allowed to change along an internal edge leading to a specified clade (see argument "clade" and Slater and Pennell in press for an example)
}

\item{eb.type }{
The type of exponential change model assumed. If eb.type = "exponential" (the default), then an exponentially declining rate will be assumed and contrasts will be log transformed when computing the node height test. If eb.type = "linear", a linear decline in rate will be assumed and untransformed contrasts will be used.
}

\item{clade }{
Default = NULL and is used if model = "BM" or model = "EB". If using model = "edge.shift", then a clade must be specified for which the stem lineage experiences a different rate of evolution. The clade is specified by giving the names of two taxa spanning the clade of interest, e.g. clade = c("A", "B") 
}

\item{rlm.maxit}{
Maximum number of interations to use for the iteratively reweighted least squares optimization of the robust regression algorithm (see ?rlm). Default is 20 and should be sufficient for most problems
}

}
\details{
This function runs a posterior predictive MCMC under the specified model, sampling model parameters from their posterior distributions and simulating under that model. Simulated data are summarized using the Node height test (Freckleton and Harvey 2006) slope  (OLS and robust regression) and Morphological Disparity Index (Harmon et al. 2003). Model adequacy can then be assessed by comparing observed values for these summary statistics to the posterior predictive distributions
}
\value{

A dataframe containing the following columns:
 \item{$generation }{the generation at which parameters where sampled and simulations conducted}
 \item{$logLk }{The sampled logLikelihood values for the model}
 \item{$Sigma }{Brownian rate parameter values}
 \item{$node.height.slope.lm }{posterior predictive distribution of slopes for the node height test using an ordinary least squares regression}
 \item{$node.height.slope.rlm }{posterior predictive distribution of slopes for the node height test using a robust regression}
 \item{$MDI }{posterior predictive distribution of MDI values}
}
\references{
Slater GJ and MW Pennell (2014) Robust regression and posterior predictive simulation increase power to detect early bursts of trait evolution. \emph{Systematic Biology}.

Freckleton RP and PH Harvey (2006) Detecting non-brownian evolution in adaptive radiations. \emph{PLoS Biology} 4:e373.

Harmon LJ, JA Schulte, A Larson, and JB Losos (2003). Tempo and mode of evolutionary readiations in iguanian lizards. \emph{Science} 301:961-964.
}

\author{
Graham Slater and Matthew Pennell
}

\seealso{
\code{\link[geiger]{nh.test}}, \code{\link[geiger]{dtt}}, \code{\link[geiger]{fitContinuous}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\donttest{
data(whales)

tmp <- treedata(whales$phy, whales$dat[,1])

phy <- tmp$phy
dat <- tmp$data[,1]

## compute observed statistics

nht.ols <- nh.test(phy, dat, regression.type = "lm",
log = TRUE, show.plot = FALSE)$coefficients[2,1]

nht.rlm <- nh.test(phy, dat, regression.type = "rlm",
 log = TRUE, show.plot = FALSE)$coefficients[2,1]

mdi.exp <- 0

#---- run short pp.mcmc 

pp.eb <- pp.mcmc(phy, dat, Ngens = 1000, sampleFreq = 10, printFreq = 100, model ="EB")

# ---- plot results

quartz(width = 5, height = 7)
par(mar = c(4,5,1,1))
par(mfcol = c(3,1))


hist(pp.eb$MDI, col = "gray", border = "gray", main = NULL, xlab = "pp.MDI", 
ylab = "Frequency", cex.axis = 1.2)

abline(v = mdi.exp, col = "black", lwd = 3, lty = 2)

mdi.p <- length(which(pp.eb$MDI<=0))/length(pp.eb$MDI)

hist(pp.eb$node.height.slope.lm, col = "gray", border = "gray", main = NULL, xlab = "pp.nht_ols", 
ylab = "Frequency", cex.axis = 1.2)

abline(v = nht.ols, col = "black", lwd = 3, lty = 2)

node.height.ols.p <- length(which(pp.eb$node.height.slope.lm <= nht.ols)) / 
(length(pp.eb$node.height.slope.lm) +1)


hist(pp.eb$node.height.slope.rlm, col = "gray", border = "gray", main = NULL, xlab = "pp.nht_ols", 
ylab = "Frequency", cex.axis = 1.2)

abline(v = nht.rlm, col = "black", lwd = 3, lty = 2)

node.height.rr.p <- length(which(pp.eb$node.height.slope.rlm <= nht.rlm)) / 
(length(pp.eb$node.height.slope.rlm) +1)

}
}

