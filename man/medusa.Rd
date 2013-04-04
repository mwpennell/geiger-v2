\name{medusa}
\alias{medusa}
\title{
MEDUSA: modeling evolutionary diversification using stepwise AIC
}
\description{
fitting piecewise diversification models according to phylogenetic and taxonomic likelihoods 
}

\usage{
medusa(phy, richness = NULL, criterion = c("aicc", "aic"),
    partitions=NA, model = c("mixed", "bd", "yule"), 
    cut = c("both","stem","node"), 
    init = c(r=0.05, epsilon=0.5), \dots)
}

\arguments{
  \item{phy}{an ultrametric phylogenetic tree of class 'phylo'}
  \item{richness}{an optional matrix of taxonomic richnesses (see \bold{Details})}
  \item{criterion}{an information criterion to use as stopping criterion}
  \item{partitions}{the number of models to be considered regardless of fit}
  \item{model}{the flavor of piecewise models to consider (see \bold{Details})}
  \item{cut}{determines where shifts are placed on the tree (see \bold{Details})}
  \item{init}{initial conditions for net diversification and relative extinction} 
  \item{...}{additional arguments to be passed to \code{\link[geiger]{treedata}}}
}

\details{
The MEDUSA model fits increasingly complex diversification models to a dataset including \code{richness} information for sampled tips in \code{phy}. The tree must have branch lengths proportional to time. The \code{richness} object is optional, but must be given if the tree is not completely sampled. MEDUSA assumes that the entire extant diversity in the group is sampled either in \code{phy} or given by information contained within the \code{richness} object. The \code{richness} object associates species richness with lineages sampled in the tree. For instance, if a genus containing a total of 10 species is exemplied in the tree by a single tip, the total diversity of the clade must be recorded in the \code{richness} object (see \bold{Examples}). All taxa missing from the tree have to be assigned to one of the tips in the \code{richness} matrix. If the \code{richness} object is \code{NULL}, the tree is assumed to be completely sampled.  

The algorithm first fits a single diversification model to the entire dataset. A series of single breakpoints in the diversification process is then added, so that different parts of the tree evolve with different parameter values (per-lineage net diversification--\code{r} and relative extinction rates--\code{epsilon}). Initial values for these diversification parameters are given through the \code{init} argument and may need to be tailored for particular datasets. The algorithm compares all single-breakpoint models to the initial model, and retains the best breakpoint. Then all possible two-breakpoint models are compared with the best single-breakpoint model, and so on. Breakpoints may be considered at a \code{"node"}, a \code{"stem"} branch, or both (as dictated by the \code{cut} argument). Birth-death or pure-birth (Yule) processes (or a combination of these processes) may be considered by the MEDUSA algorithm. The model flavor is determined through the \code{model} argument. 

Two stopping criteria are available for the MEDUSA algorithm. The user can either determine the number of models explored by MEDUSA or this number may be determined based on model fits. A number of model partitions to be explored may be given a priori (if given a non-zero number through the \code{partitions} argument) or an information \code{criterion} is used to choose sufficient model complexity. If the latter stopping criterion is used, one needs to specify whether to use Akaike information criterion (\code{"aic"}) or sample-size corrected AIC (\code{"aicc"}). The algorithm defaults to using AICc. An appropriate threshold in AICc differences between different MEDUSA models has been shown to be dependent on tree size. The threshold used for model selection is computed internally (and is based on extensive simulation study). This value is reported to the user after the algorithm finishes. 

The user will almost certainly want to summarize the object returned from MEDUSA with the function, \code{\link[geiger]{summary.medusaRAW}}. 
}

\value{
A list object is returned including fits for all model complexities as well as summary information:

\item{control}{is a list object specifying the stopping criterion, the information criterion and threshold used (if appropriate), or the number of partitions explored}
\item{cache}{is a list object primarily used internally (including the tree and richness information)}
\item{models}{is a list object containing each optimized piecewise model, primarily for internal use}
\item{summary}{is a dataframe containing breakpoints and fit values for optimal models at each model complexity. If \code{partitions} is used as a stopping criterion, AIC and AICc weights (\code{w.aic} and \code{w.aicc}) for each model are reported (otherwise, these weights are meaningless). Other data include: the number of parameters for each model (\code{k}, determined by the number of breakpoints, independent net-diversification rates, and independent relative-extinction values), the branch or node chosen for each successive piecewise model (\code{split}), whether the split occurs at a node or stem (\code{cut}), and the model likelihood (\code{lnL})}
\item{FUN}{is a function used to summarize a particular model (indexed by number) and is primarily for internal use}
}


\references{
Alfaro, ME, F Santini, C Brock, H Alamillo, A Dornburg, DL Rabosky, G Carnevale, and LJ Harmon. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. \emph{Proceedings of the National Academy of Sciences} \bold{106}: 13410-13414.
}

\author{
JW Brown <phylo.jwb@gmail.com>, RG FitzJohn, ME Alfaro, LJ Harmon, and JM Eastman
}

\seealso{
\code{\link[geiger]{summary.medusaRAW}}
\code{\link[geiger]{plot.medusaRAW}}
\code{\link[geiger]{plot.medusa}}
}

\examples{
    \donttest{
    dat=get(data(whales))
    phy=dat$phy
    richness=dat$richness
    
    ## USING AICc as STOPPING CRITERION
    res1=medusa(phy, richness, warnings=FALSE)
    print(names(res1)) # output list elements
    print(res1$summary) # show 'summary' object
    summary(res1, criterion="aicc") # select best model based on AICc
    
    ## USING PARTITIONS as STOPPING CRITERION

    res2=medusa(phy, richness, partitions=12, warnings=FALSE)
    s2=summary(res2, criterion="aicc") # select best model based on AICc (again showing the third model as best)
    
    ## PLOTTING RESULTS
    # plot model fit by model complexity for AIC and AICc
    plot(res2, col=c(aic="gray", aicc="black")) # invoking plot.medusaRAW()
    
    # plot breakpoints for the best model chosen by AICc
    plot(s2, cex=0.5, partitions=list(cex=5), label.offset=1, edge.width=2) # invoking plot.medusa()
    }
}
\keyword{models}


