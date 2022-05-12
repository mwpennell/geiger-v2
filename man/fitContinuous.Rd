\name{fitContinuous}
\alias{fitContinuous}
\title{ Model fitting for continuous comparative data }
\description{
	fitting macroevolutionary models to phylogenetic trees
}
\usage{
fitContinuous(phy, dat, SE = 0,
    model = c("BM","OU","EB","rate_trend","lambda","kappa","delta","mean_trend","white"),
    bounds= list(), control = list(method = c("subplex","L-BFGS-B"),
    niter = 100, FAIL = 1e+200, hessian = FALSE, CI = 0.95), ncores=NULL, \dots)
}
\arguments{
  \item{phy}{ a phylogenetic tree of class phylo}
  \item{dat}{ data vector for a single trait, with names matching tips in \code{phy}}
  \item{SE}{ a single value or named vector of standard errors associated with values in \code{dat}; if any elements in the vector \code{SE} are NA, SE will be estimated }
  \item{model}{ model to fit to comparative data (see \bold{Details}) }
  \item{bounds}{ range to constrain parameter estimates (see \bold{Details})  }
  \item{control}{ settings used for optimization of the model likelihood}
  \item{ncores}{ Number of cores. If \code{NULL} then number of cores is detected}
  \item{...}{additional arguments to be passed to the internal likelihood function \code{bm.lik}}
}
\details{
This function fits various likelihood models for continuous character evolution. The function
returns parameter estimates and the likelihood for univariate datasets.

The model likelihood is maximized using methods available in \code{\link[stats]{optim}} as well as \code{\link[subplex]{subplex}}.
Optimization methods to be used within \code{optim} can be specified through the \code{control} object (i.e., \code{control$method}).

A number of random starting points are used in optimization and are given through the \code{niter} element within the \code{control} object
(e.g., \code{control$niter}). The \code{FAIL} value within the \code{control} object should be a large value
that will be considerably far from -lnL of the maximum model likelihood. In most cases, the default setting
for \code{control$FAIL} will be appropriate. The Hessian may be used to compute confidence intervals
(\code{CI}) for the parameter estimates if the \code{hessian} element in \code{control} is TRUE.

\emph{Beware}: difficulty in finding the optimal solution is determined by an interaction between the nature and complexity of the
likelihood space (which is data- and model-dependent) and the numerical optimizer used to explore the space.  There is never a
guarantee that the optimal solution is found, but using many random starting points (\code{control$niter}) and many optimization methods
(\code{control$method}) will increase these odds.

Bounds for the relevant parameters of the fitted model may be given through the \code{bounds} argument. Bounds
may be necessary (particularly under the \code{OU} model) if the likelihood surface is characterized by a long,
flat ridge which can be exceedingly difficult for optimization methods. Several bounds can be given at a time
(e.g., \code{bounds=list(SE=c(0,0.1),alpha=c(0,1))} would constrain measurement error as well as the 'constraint'
parameter of the Ornstein-Uhlenbeck model).  Default bounds under the different models are given below.

Possible models are as follows:
\itemize{
	\item{\bold{BM} }{is the Brownian motion model (Felsenstein 1973), which assumes the correlation structure among trait values is proportional to the extent of shared ancestry
		for pairs of species. Default bounds on the rate parameter are \code{sigsq=c(min=exp(-500),max=exp(100))}. The same bounds are applied to all other
		models, which also estimate \code{sigsq}}
	\item{\bold{OU} }{is the Ornstein-Uhlenbeck model (Butler and King 2004), which fits a random walk with a central tendency with an attraction strength proportional to the parameter \code{alpha}.
		The \code{OU} model is called the \code{hansen} model in \pkg{ouch}, although the way the parameters are fit is slightly different here. Default
		bounds are \code{alpha = c(min = exp(-500), max = exp(1))}}
	\item{\bold{EB} }{is the Early-burst model (Harmon et al. 2010) and also called the \code{ACDC} model (accelerating-decelerating; Blomberg et al. 2003). Set by the \code{a} rate parameter, \code{EB}
		fits a model where the rate of evolution increases or decreases exponentially through time, under the model r[t] = r[0] * exp(a * t), where \code{r[0]} is the
		initial rate, \code{a} is the rate change parameter, and \code{t} is time. The maximum bound is set to \code{-0.000001}, representing a decelerating rate of evolution. The minimum bound is set to log(10^-5)/depth of the tree.}
	\item{\bold{rate_trend} }{is a diffusion model with linear trend in rates through time (toward larger or smaller rates). Used to be denominated the \code{"trend"} model, which is still accepted by \code{fitContinuous} for backward compatibility. Default bounds are \code{slope = c(min = -100, max = 100)}}
	\item{\bold{lambda} }{is one of the Pagel (1999) models that fits the extent to which the phylogeny predicts covariance among trait values for species. The model effectively transforms the tree:
		values of \code{lambda} near 0 cause the phylogeny to become more star-like, and a \code{lambda} value of 1 recovers the \code{BM} model. Default
		bounds are \code{lambda = c(min = exp(-500), max = 1}}
	\item{\bold{kappa} }{is a punctuational (speciational) model of trait evolution (Pagel 1999), where character divergence is related to the number of speciation events between two species. Note that if
		there are speciation events that are missing from the given phylogeny (due to extinction or incomplete sampling), interpretation under the \code{kappa} model may be difficult. Considered as a tree
		transformation, the model raises all branch lengths to an estimated power (\code{kappa}). Default bounds are \code{kappa = c(min = exp(-500), max = 1)}  }
	\item{\bold{delta} }{is a time-dependent model of trait evolution (Pagel 1999).  The \code{delta} model is similar to \code{ACDC} insofar as the \code{delta} model fits the relative contributions of
		early versus late evolution in the tree to the covariance of species trait values.  Where \code{delta} is greater than 1, recent evolution has been relatively fast; if \code{delta} is less
		than 1, recent evolution has been comparatively slow. Intrepreted as a tree transformation, the model raises all node depths to an estimated power (\code{delta}). Default bounds are \code{delta = c(min = exp(-500), max = 3)}}
    \item{\bold{mean_trend} }{is a model of trait evolution with a directional drift or \emph{trend} component (i.e., toward smaller or larger values through time).  This model is sensible only for non-ultrametric trees, as the likelihood surface is entirely flat with respect to the slope of the trend if the tree is ultrametric. The model used to be denominated the \code{"drift"} model, which is still accepted by \code{fitContinuous} for backward compatibility. Default bounds are \code{drift = c(min = -100, max = 100)}}
	\item{\bold{white} }{is a \code{white}-noise (non-phylogenetic) model, which assumes data come from a single normal distribution with no covariance structure among species. The variance parameter \code{sigsq}
		takes the same bounds defined under the \code{BM} model}
}
}
\value{
	\code{fitContinuous} returns a list with the following four elements:

		\item{\bold{lik} }{is the function used to compute the model likelihood. The returned function (\code{lik}) takes arguments that are necessary for the given model.
			For instance, if estimating a Brownian-motion model with unknown standard error, the arguments (\code{pars}) to the \code{lik} function would be \code{sigsq} and \code{SE}. By default, the function evaluates the likelihood of the model by assuming the maximum likelihood root state. This behavior can be changed in the call to \code{lik} with \code{lik(pars, root=ROOT.GIVEN)} where \code{pars} includes a value for the root state (z0). See \bold{Examples} for a demonstration. The tree and data are stored internally within the \code{lik} function, which permits those elements to be efficiently reused when computing the likelihood under different parameter values
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
\note{
    To speed the likelihood search, one may set an environment variable to make use of parallel processing, used by \code{\link[parallel]{mclapply}}. To set the environment variable, use \code{options(mc.cores=INTEGER)}, where \code{INTEGER} is the number of available cores. Alternatively, the \code{mc.cores} variable may be preset upon the initiation of an R session (see \code{\link[base]{Startup}} for details).
}
\references{
Blomberg SP, T Garland, and AR Ives. 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. \emph{Evolution} 57:717-745.

Butler MA and AA King, 2004. Phylogenetic comparative analysis: a modeling approach for adaptive evolution. \emph{American Naturalist} 164:683-695.

Felsenstein J. 1973. Maximum likelihood estimation of evolutionary trees from continuous characters. \emph{American Journal of Human Genetics} 25:471-492.

Harmon LJ et al. 2010. Early bursts of body size and shape evolution are rare in comparative data. \emph{Evolution} 64:2385-2396.

Pagel M. 1999. Inferring the historical patterns of biological evolution. \emph{Nature} 401:877-884
}

\author{LJ Harmon, W Challenger, and JM Eastman}
\examples{
\dontrun{
geo=get(data(geospiza))

tmp=treedata(geo$phy, geo$dat)
phy=tmp$phy
dat=tmp$data

#---- STORE RESULTS
brownFit <-  fitContinuous(phy, dat[,"wingL"], SE=NA, control=list(niter=50), ncores=2)

#---- PRINT RESULTS
print(names(brownFit))
print(brownFit)

\donttest{
#---- COMPUTE LIKELIHOOD
flik=brownFit$lik
print(argn(flik))

#---- CREATE a FUNCTION to COMPARE MODELS
fitGeospiza=function(trait=c("wingL","tarsusL","culmenL","beakD","gonysW")){

	trait=match.arg(trait, c("wingL","tarsusL","culmenL","beakD","gonysW"))

	# define set of models to compare
	models=c("BM", "OU", "EB", "white")
	summaries=c("diffusion", "Ornstein-Uhlenbeck", "early burst", "white noise")

	## ESTIMATING measurement error ##
	aic.se=numeric(length(models))
	lnl.se=numeric(length(models))

	for(m in 1:length(models)){
		cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m],
			" with SE *** \n", sep="")
		tmp=fitContinuous(phy,dat[,trait],SE=NA, model=models[m],
                                    bounds=list(SE=c(0,0.5)), ncores=2)
		print(tmp)
		aic.se[m]=tmp$opt$aicc
		lnl.se[m]=tmp$opt$lnL
	}


	## ASSUMING no measurement error ##
	aic=numeric(length(models))
	lnl=numeric(length(models))

	for(m in 1:length(models)){
		cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m],
			 " *** \n", sep="")
		tmp=fitContinuous(phy,dat[,trait],SE=0,model=models[m], ncores=2)
		print(tmp)
		aic[m]=tmp$opt$aicc
		lnl[m]=tmp$opt$lnL
	}

	## COMPARE AIC ##
	names(aic.se)<-names(lnl.se)<-names(aic)<-names(lnl)<-models
	delta_aic<-function(x) x-x[which(x==min(x))]

	# no measurement error
	daic=delta_aic(aic)
	cat("\n\n\n\t\t\t\t*** MODEL COMPARISON: ",trait," *** \n",sep="")
	cat("\tdelta-AIC values for models assuming no measurement error
    \t\t\t\t zero indicates the best model\n\n")
	print(daic, digits=2)

		# measurement error
	daic.se=delta_aic(aic.se)
	cat("\n\n\n\n\t\t\t\t*** MODEL COMPARISON: ",trait," ***\n",sep="")
	cat("\t\t   delta-AIC values for models estimating SE
    \t\t\t\t zero indicates the best model\n\n")
	print(daic.se, digits=2)
	cat("\n\n\n")

	res_aicc=rbind(aic, aic.se, daic, daic.se)
	rownames(res_aicc)=c("AICc","AICc_SE","dAICc", "dAICc_SE")

	return(res_aicc)
}

#---- COMPARE MODELS for WING LENGTH
res=fitGeospiza("wingL")
print(res)
}
}
}

\keyword{models}
