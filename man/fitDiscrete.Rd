\name{fitDiscrete}
\alias{fitDiscrete}
\title{ model fitting for discrete comparative data }
\description{
	fitting macroevolutionary models to phylogenetic trees
}
\usage{
fitDiscrete(phy, dat,
    model = c("ER","SYM","ARD","meristic"),
    transform = c("none", "EB", "lambda", "kappa", "delta", "white"),
    bounds = list(), control = list(method = c("SANN", "L-BFGS-B"),
    niter = 100, FAIL = 1e+200, hessian = FALSE, CI = 0.95), ...)
}
\arguments{
  \item{phy}{ a phylogenetic tree of class phylo}
  \item{dat}{ data vector for a single trait, with names matching tips in \code{phy}}
  \item{model}{ an Mkn model to fit to comparative data (see \bold{Details}) }
  \item{transform}{ an evolutionary model used to transform the tree (see \bold{Details}) }
  \item{bounds}{ range to constrain parameter estimates (see \bold{Details})  }
  \item{control}{ settings used for optimization of the model likelihood}
  \item{...}{if \code{model="meristic"}, \code{...} can dictate whether the matrix is asymmetric (\code{symmetric=FALSE})}
}
\details{
This function fits various likelihood models for discrete character evolution. The function 
returns parameter estimates and the likelihood for univariate datasets. 

The model likelihood is maximized using methods available in \code{\link[stats]{optim}} as well as \code{\link[subplex]{subplex}}. 
Optimization methods to be used within \code{optim} can be specified through the \code{control} object. A number of random starting points 
are used in optimization and are given through the \code{niter} element within the \code{control} object 
(e.g., \code{control$niter}). The \code{FAIL} value within the \code{control} object should be a large value 
that will be considerably far from -lnL of the maximum model likelihood. In most cases, the default setting 
for \code{control$FAIL} will be appropriate. The Hessian may be used to compute confidence intervals 
(\code{hessian_CI}) for the parameter estimates if the \code{hessian} element in \code{control} is TRUE.

All of the models are continuous-time Markov models of trait evolution (see Yang 2006 for a good general discussion of this type of model).  
The function can handle traits with any number of character states, under a range of models. The character model is specified by the \code{model} argument:

\itemize{
	\item{\bold{ER} }{is an \code{equal-rates} model of where a single parameter governs all transition rates}
	\item{\bold{SYM} }{is a \code{symmetric} model where forward and reverse transitions share the same parameter}
	\item{\bold{ARD} }{is an \code{all-rates-different} model where each rate is a unique parameter}  
	\item{\bold{meristic} }{is a model wherein transitions occur in a stepwise fashion (e.g., 1 to 2 to 3 to 2) without skipping intermediate steps; this requires a sensible coding of the character 
		states as consecutive integers are assumed to be neighboring states}
	\item{\bold{matrix} }{is a user supplied model (given as a dummy matrix representing transition classes between states); elements that are zero signify rates that are also zero (see \bold{Examples})}
}


The \code{transform} argument allows one to test models where rates vary across the tree.  Bounds for the relevant parameters of the tree \code{transform} 
may be given through the \code{bounds} argument. Several bounds can be given at a time. Default bounds under the different models are given below. 
Options for \code{transform} are as follows:

\itemize{
	\item{\bold{none} }{is a model of rate constancy through time}
	\item{\bold{EB} }{is the Early-burst model (Harmon et al. 2010) and also called the \code{ACDC} model (accelerating-decelerating; Blomberg et al. 2003). Set by the \code{a} rate parameter, \code{EB} 
		fits a model where the rate of evolution increases or decreases exponentially through time, under the model r[t] = r[0] * exp(a * t), where \code{r[0]} is the 
		initial rate, \code{a} is the rate change parameter, and \code{t} is time. Default bounds are \code{a = c(min = -10, max = 10)}}
	\item{\bold{lambda} }{is one of the Pagel (1999) models that fits the extent to which the phylogeny predicts covariance among trait values for species. The model effectively transforms the tree: 
		values of \code{lambda} near 0 cause the phylogeny to become more star-like, and a \code{lambda} value of 1 recovers the \code{none} model. Default 
		bounds are \code{lambda = c(min = 0, max = 1)}}
	\item{\bold{kappa} }{is a punctuational (speciational) model of trait evolution (Pagel 1999), where character divergence is related to the number of speciation events between two species. Note that if 
		there are speciation events in the given phylogeny (due to extinction or incomplete sampling), interpretation under the \code{kappa} model may be difficult. Considered as a tree 
		transformation, the model raises all branch lengths to an estimated power (\code{kappa}). Default bounds are \code{kappa = c(min = 0, max = 1)}}
	\item{\bold{delta} }{is a time-dependent model of trait evolution (Pagel 1999).  The \code{delta} model is similar to \code{ACDC} insofar as the \code{delta} model fits the relative contributions of 
		early versus late evolution in the tree to the covariance of species trait values.  Where \code{delta} is greater than 1, recent evolution has been relatively fast; if \code{delta} is less 
		than 1, recent evolution has been comparatively slow. Intrepreted as a tree transformation, the model raises all node depths to an estimated power (\code{delta}). Default bounds are \code{delta = c(min = 0, max = 3)}}
	\item{\bold{white} }{is a \code{white}-noise (non-phylogenetic) model, which converts the tree into a star phylogeny}
}
}

\value{
	\code{fitDiscrete} returns a list with the following four elements: 
	
		\item{\bold{lik} }{is the function used to compute the model likelihood. The returned function (\code{lik}) takes arguments that are necessary for the given model. 
			For instance, if estimating an untransformed \code{ER} model, there would be a single argument (the transition rate) necessary for the \code{lik} function. The tree and data are stored internally within the \code{lik} function, which permits those elements to be efficiently reused when computing 
			the likelihood under different parameter values 
		}
		\item{\bold{bnd} }{is a matrix of the used bounds for the relevant parameters estimated in the model. Warnings will be issued if any parameter estimates occur at the supplied (or default) parameter bounds
		}
		\item{\bold{res} }{is a matrix of results from optimization. Rownames of the \code{res} matrix are the optimization methods 
			(see \code{\link[stats]{optim}} and \code{\link[subplex]{subplex}}). The columns in the \code{res} matrix are the estimated 
			parameter values, the estimated model likelihood, and an indication of optimization convergence. Values of convergence not 
			equal to zero are not to be trusted  
		}
		\item{\bold{opt} }{is a list of the primary results: estimates of the parameters, the maximum-likelihood estimate (\code{lnL}) of the model, the 
			optimization method used to compute the MLE, the number of model parameters (\code{k}, including one parameter for the root state), the AIC (\code{aic}), 
			sample-size corrected AIC (\code{aicc}). The number of observations for AIC computation is taken to be the number of trait values observed. 
			If the Hessian is used, confidence intervals on the parameter estimates (\code{CI}) and the Hessian matrix (\code{hessian}) are also returned 
		}
}
\references{ 
Yang, Z. 2006. \emph{Computational Molecular Evolution}. Oxford University Press: Oxford.
}
\author{ LJ Harmon, RE Glor, RG FitzJohn, and JM Eastman }
\examples{


## match data and tree
tmp=get(data(geospiza))
td=treedata(tmp$phy, tmp$dat)
geo=list(phy=td$phy, dat=td$data)
gb=round(geo$dat[,5]) ## create discrete data
names(gb)=rownames(geo$dat)
mm=rbind(c(NA, 0, 0), c(1, NA, 2), c(0, 2, NA))

fitDiscrete(geo$phy, gb, model="ER", control=list(niter=5)) #-7.119792
\donttest{
fitDiscrete(geo$phy, gb, model="ARD") #-6.064573
fitDiscrete(geo$phy, gb, model=mm) #-7.037944
fitDiscrete(geo$phy, gb, model="SYM", transform="kappa")} #-6.822943

}
\keyword{ arith }
\keyword{models}