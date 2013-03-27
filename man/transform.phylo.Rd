\name{transform.phylo}
\alias{transform.phylo}

\title{Tree transformations}
\description{
	applying various transformation to the branches of a phylogenetic tree.
	}
\usage{
\method{transform}{phylo}(x, model = c("OU", "EB", "nrate", "lrate", 
    "trend", "lambda", "kappa", "delta", "white", "depth"), \dots)
}
\arguments{
  \item{x}{an object of class phylo }
  \item{model}{a model used to transform the tree (see \bold{Details})}
  \item{...}{argument(s) to be passed to the transformation function (see \bold{Examples})}
}
\details{
	This function takes a tree and returns either a transformed tree if \code{...} is not empty and gives the parameter value(s) for the tree 
	transformation. If \code{...} is left empty, a function is returned to the user that can be efficiently iterated over many parameter values 
	for transformation. The available models are meant to correspond with changing the model of phenotypic evolution 
	for discrete or continuous characters. 
}
\value{
	A transformation function (or transformed phylogenetic tree of class 'phylo' (\pkg{ape} format) is returned. Possible transforms include the following: 
    
        \item{\bold{OU} }{is the Ornstein-Uhlenbeck model (Butler and King 2004), which fits a random walk with a central tendency with an attraction strength proportional to the parameter \code{alpha}.
		The \code{OU} model is called the \code{hansen} model in \pkg{ouch}, although the way the parameters are fit is slightly different here. The parameter used for transformation is \code{alpha}}

		\item{\bold{EB} }{is the Early-burst model (Harmon et al. 2010) and also called the \code{ACDC} model (accelerating-decelerating; Blomberg et al. 2003). Set by the \code{a} rate parameter, \code{EB} 
		fits a model where the rate of evolution increases or decreases exponentially through time, under the model r[t] = r[0] * exp(a * t), where \code{r[0]} is the 
		initial rate, \code{a} is the rate change parameter, and \code{t} is time. The parameter used for transformation is \code{a}.}

        \item{\bold{nrate} }{is the multiple-rates model where time slices have independent rates of evolution. The parameters used for transformation are \code{time} and \code{rate}, both of which may be vectors. If rates for several time slices are given, the vectors \code{time} and \code{rate} must match in length. The \code{time} argument is expected in relative branching heights (where the root is 0 and the tips of an ultrametric tree terminate at a relative time of 1). Note that the default behavior is for the multirate transformation to rescale the tree to its original height, while preserving the relative rates across the tree. Note also that the initial relative rate (for the rootmost time slice) is 1.}

        \item{\bold{lrate} }{is the multiple-rates model where local clades have independent rates of evolution. The parameters used for transformation are \code{node} and \code{rate}, both of which may be vectors. If rates for several lineages are given, the vectors \code{node} and \code{rate} must match in length. The \code{node} argument is expected to have node identifiers consistent with \pkg{ape} labeling (see \code{\link[ape]{read.tree}}). This \code{node} argument defines where breakpoints occur in the tree (at which a transition to the associated relative rate occurs). Note that the stem branch associated with the node is included as part of the defined rate partition. Rates given are relative to an assumed rootmost partition with a rate scalar of 1.}

		\item{\bold{trend} }{is a diffusion model with linear trend in character value through time (toward larger or smaller trait values). Support for the \code{trend} model will be lacking for 
		ultrametric trees. The parameter used for transformation is \code{slope}.}

		\item{\bold{lambda} }{is one of the Pagel (1999) models that fits the extent to which the phylogeny predicts covariance among trait values for species. The model effectively transforms the tree: 
		values of \code{lambda} near 0 cause the phylogeny to become more star-like, and a \code{lambda} value of 1 recovers the \code{BM} model. The parameter used for transformation is \code{lambda}.}

		\item{\bold{kappa} }{is a punctuational (speciational) model of trait evolution (Pagel 1999), where character divergence is related to the number of speciation events between two species. Note that if 
		there are missing speciation events in the given phylogeny (due to extinction or incomplete sampling), interpretation under the \code{kappa} model may be difficult. Considered as a tree 
		transformation, the model raises all branch lengths to an estimated power (\code{kappa}). The parameter used for transformation is \code{kappa}.}

		\item{\bold{delta} }{is a time-dependent model of trait evolution (Pagel 1999).  The \code{delta} model is similar to \code{ACDC} insofar as the \code{delta} model fits the relative contributions of 
		early versus late evolution in the tree to the covariance of species trait values.  Where \code{delta} is greater than 1, recent evolution has been relatively fast; if \code{delta} is less 
		than 1, recent evolution has been comparatively slow. Intrepreted as a tree transformation, the model raises all node depths to an estimated power (\code{delta}). The parameter used for transformation 
		is \code{delta}. Note that the default behavior is for the delta transformation to rescale the tree to its original height.}   

		\item{\bold{white} }{is a \code{white}-noise (non-phylogenetic) model, which assumes data come from a single normal distribution with no covariance structure among species}

		\item{\bold{depth} }{is simply a transformation of the total depth of the tree; stretching the tree has an effect of increasing rates of evolution under Brownian motion (relative to characters evolved 
		on the unstretched tree), and compressing the tree has the opposite effect. The parameter used for transformation is \code{depth}.}

}
\references{
Pagel, M. 1999. Inferring the historical patterns of biological evolution. \emph{Nature} 401:877-884.
			
Butler, M.A. and A.A. King, 2004. Phylogenetic comparative analysis: a modeling approach for adaptive evolution. \emph{American Naturalist} 164:683-695.

Various papers in prep., L. J. Harmon and J. T. Weir.

}
\author{LJ Harmon and JM Eastman}
\examples{


geo=get(data(geospiza))

## returning a function
ltrns=transform(geo$phy, "lambda")
plot(ltrns(0))
title("lambda: 0.0")

plot(ltrns(0.5))
title("lambda: 0.5")

plot(ltrns(1))
title("lambda: 1")


## transforming the tree
lphy=transform(geo$phy, "lambda", 0.5) # transform tree in one fell swoop
plot(lphy)
title("lambda: 0.5")

## multirate tree -- time
rtrns=transform(geo$phy, "nrate")
rphy=rtrns(time=c(0.2, 0.4, 0.6, 0.8), rate=c(2, 4, 8, 16), rescale=TRUE)
plot(rphy)
title("5-rate tree: by time")

## multirate tree -- lineages
mtrns=transform(geo$phy, "lrate")
mphy=mtrns(node=c(25, 20), rate=c(4, 8))
plot(mphy)
title("3-rate tree: by lineages")

}
\keyword{manip}
\keyword{ graphs }
