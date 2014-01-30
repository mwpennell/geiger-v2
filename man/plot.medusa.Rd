\name{plot.medusa}
\alias{plot.medusa}
\title{
MEDUSA: modeling evolutionary diversification using stepwise AIC
}
\description{
summarizing piecewise diversification models estimated by MEDUSA 
}

\usage{
	\method{plot}{medusa}(x, cex = 0.5, time = TRUE, ...)
}

\arguments{
  \item{x}{an object of class \code{medusa} for plotting}
  \item{cex}{text size}
  \item{time}{logical. should a time axis be plotted. default = TRUE.}
  \item{...}{additional arguments to be passed to internal functions}
}

\details{
The \code{\link[geiger]{medusa}} model returns a raw list object. This function is used to generate a modified \code{edge} matrix (see \code{\link[ape]{read.tree}} for details on the \code{edge} matrix), giving all relevant information about the estimated diversification process. The returned \bold{z}-matrix includes: the ancestor (\code{anc}) and descendant (\code{dec}) relationships between nodes of the tree (using \pkg{ape} indices); the beginning (\code{t.0}) and ending (\code{t.1}) times and length (\code{t.len}) of each branch; the diversities at the start (\code{n.0}) and end (\code{n.t}) of each branch; the piecewise model assigned \code{partition} to the branch; whether the branch is associated with a \code{shift}; the timing of the shift (\code{t.shift}); the net-diversification rate (\code{r}) and relative-extinction rate (\code{epsilon}) associated with the branch as well as for the direct ancestor of the branch (\code{ancestral.r} and \code{ancestral.epsilon}). The  \bold{z}-matrix also includes a \code{summary} attribute that shows which model is chosen and associated information on model fit (see \bold{Examples}). 

The raw output of \code{\link[geiger]{medusa}} contains an optimized model. The summary output may then be sent to a plotting function which will display the location on the tree where breakpoints have been placed. Note that the first piecewise model corresponds to the root and all descendants (until another breakpoint is encountered).
}

\references{
Alfaro, ME, F Santini, C Brock, H Alamillo, A Dornburg, DL Rabosky, G Carnevale, and LJ Harmon. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. \emph{Proceedings of the National Academy of Sciences} \bold{106}: 13410-13414.
}

\author{
JW Brown <phylo.jwb@gmail.com>, RG FitzJohn, ME Alfaro, LJ Harmon, and JM Eastman
}

\seealso{
\code{\link[geiger]{medusa}}
}

\examples{
    \donttest{
    dat=get(data(whales))
    phy=dat$phy
    richness=dat$richness
        
    res <- medusa(phy, richness)
    
    # select best model based on AICc (showing the third model as best)
    plot(res, cex=0.5, label.offset=1) # using plot.medusa()
    title("AICc-chosen model")
    }
}
\keyword{graphics}

