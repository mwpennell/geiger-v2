# optimization for single-rate Brownian model -- use fitContinuous instead 
.TESTING.likfx.bm=function(phy, dat, SE=NULL, root=NA, method=c("direct","vcv","reml")){
	method=match.arg(method, c("direct","vcv","reml"))
	
	lik=make.bm.relaxed(phy, dat, SE, method=method)
	
	if(!is.na(root) || method=="reml") {
		x=function(par) lik(rep(exp(par),nrow(attributes(lik)$cache$phy$edge)), root)
		par=c(s2=-1)
		method="Brent"
		lower=-100
		upper=100
	} else {
		x=function(par) lik(rep(exp(par[1]),nrow(attributes(lik)$cache$phy$edge)), par[2])
		par=c(s2=-1,root=0)
		method="Nelder-Mead"
		lower=-Inf
		upper=Inf
	}
	z=optim(par, fn=x, lower=lower, upper=upper, method=method, control=list(fnscale=-1))
	z$par[1]=exp(z$par[1])
	z
}

## replacing prepare.data.bm
make.bm.relaxed <- function(phy, dat, SE=NA, method=c("direct","vcv","reml")){
	method=match.arg(method, c("direct","vcv","reml"))
	lik=switch(method,
			   direct=.make.bm.relaxed.direct(phy, dat, SE),
			   vcv=.make.bm.relaxed.vcv(phy, dat, SE),
			   reml=.make.bm.relaxed.pic(phy, dat, SE))
	lik
}


#primary function for computing the likelihood of data, given a root state, VCV matrix, and Brownian motion model 
#author: LJ HARMON 2009 and JM EASTMAN 2010
.bm.lik.fn.vcv <-
function(root, dat, vcv, SE) { 
# mod 12.02.2010 JM Eastman: using determinant()$modulus rather than det() to stay in log-space
	y=dat
	b <- vcv
	if(any(SE!=0)) diag(b)=diag(b)+SE^2
	w <- rep(root, nrow(b))
	num <- -t(y-w)%*%solve(b)%*%(y-w)/2
	den <- 0.5*(length(y)*log(2*pi) + as.numeric(determinant(b)$modulus))
	return((num-den)[1,1])
}


#primary function for computing the likelihood of data, using REML, under Brownian motion model 
#author: LJ HARMON, LJ REVELL, and JM EASTMAN 2011
#'phy' and 'rates' must be in 'pruningwise' order (as 'ic')
.bm.lik.fn.reml <- 
function(phy, rates, ic) {
	new.tre=.scale.brlen(phy, rates)
	new.var=.pic_variance.phylo(new.tre)
	reml=dnorm(ic, mean=0, sd=sqrt(new.var), log=TRUE)
	return(sum(reml))	
}


#compute expected PIC variance given tree: used for .bm.lik.fn.reml()  
#author: JM EASTMAN 2011
.pic_variance.phylo <- function(phy)
{
	phy$node.label=NULL
	nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
	if(!is.binary.tree(phy)) stop("'phy' is not fully dichotomous.")
	
    phy <- reorder(phy, "pruningwise")
	
	ans <- .C("pic_variance", as.integer(nb.tip), as.integer(nb.node),
              as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
              as.double(phy$edge.length), double(nb.node),
              PACKAGE = "geiger")
	
    var <- ans[[6]]
    var
}

.make.bm.relaxed.vcv <- function(phy, dat, SE=NULL){
	cache=.prepare.bm(phy, dat, SE)
	
	check.argn=function(rates, root){
		if(length(rates)!=(cache$n.node+cache$n.tip-1) && length(root)==1) stop("Supply 'rates' as a vector of rate scalars for each branch, and supply 'root' as a single value.")
	}
	
	lik <- function(rates, root){
		check.argn(rates, root)
		tt=.scale.brlen(cache$phy, rates)
		vcv=.vmat(tt)
		SE=SE[match(colnames(vcv),names(cache$SE))]
		dat=dat[match(colnames(vcv),names(cache$dat))]
		.bm.lik.fn.vcv(root, cache$dat, vcv, cache$SE)
	}
	attr(lik,"cache")=cache
	attr(lik,"argn")=list(rates=cache$phy$edge[,2], root="root")
	class(lik)=c("rbm","bm","function")
	lik
}

.make.bm.relaxed.pic <- function(phy, dat, SE=NULL){
	cache=.prepare.bm(phy, dat, SE)
	
	ic=pic(dat, phy, scaled=FALSE)
	cache$ic=ic
	
	check.argn=function(rates, root){
		if(length(rates)!=(cache$n.node+cache$n.tip-1)) stop("Supply 'rates' as a vector of rate scalars for each branch.")
	}
	
	lik <- function(rates, root=NA){
		check.argn(rates, root)
		mm=match(cache$phy$tip.label, names(cache$SE))
		tt=cache$phy$edge[,2]<=cache$n.tip
		rates[tt]=rates[tt]+cache$SE[mm]^2
		.bm.lik.fn.reml(cache$phy, rates, cache$ic)
	}
	attr(lik,"cache")=cache
	attr(lik,"argn")=list(rates=cache$phy$edge[,2], root=NA)
	class(lik)=c("rbm","bm","function")
	lik
}

.make.bm.relaxed.direct <- function (phy, dat, SE=NULL) 
{
	cache=.prepare.bm(phy, dat, SE)
	z = length(cache$len)
    rr = numeric(z)
    rootidx = as.integer(cache$root)
	adjvar = as.integer(cache$y$SE)
	
    datc_init = list(len = as.numeric(cache$len),
				intorder = as.integer(cache$order[-length(cache$order)]), 
				tiporder = as.integer(cache$y$target), 
				root = rootidx, 
				y = as.numeric(cache$y$y[1, ]), 
				var = as.numeric(cache$y$y[2, ]),
				n = as.integer(z), 
				descRight = as.integer(cache$children[, 1]),
				descLeft = as.integer(cache$children[, 2]),
#               thetasq=0,
#               sigsq=0,
#               alpha=0,
#               hsq=0,
                model=0)
	
	ll.bm.direct <- function(pars, root = ROOT.MAX, root.x = NA, intermediates = FALSE, datc) {
        out = .Call("bm_direct", dat = datc, pars = pars, package = "geiger")
        vals = c(out$initM[rootidx], out$initV[rootidx], out$lq[rootidx])
        loglik <- .root.bm.direct(vals, out$lq[-rootidx], root, root.x)
        if (intermediates) {
            attr(loglik, "intermediates") <- intermediates
            attr(loglik, "vals") <- vals
        }
        return(loglik)
    }
    class(ll.bm.direct) <- c("bm.direct", "bm", "function")
	
	N = cache$n.tip
    n = cache$n.node
    ord = 1:(N + n)
	vv = numeric(length(ord))
    mm = match(cache$phy$edge[, 2], ord)	
	check.argn=function(rates, root){
		if(length(rates)!=(cache$n.node+cache$n.tip-1) && length(root)==1) stop("Supply 'rates' as a vector of rate scalars for each branch, and supply 'root' as a single value.")
	}
	
	## LIKELIHOOD FUNCTION
	if(any(adjvar==1)){ # adjustable SE
		likSE <- function(rates, root, SE) {
			check.argn(rates, root)
			vv[mm] = rates
			datc_se=datc_init
			datc_se$var[which(adjvar==1)]=SE^2
			ll = ll.bm.direct(pars = vv, root = ROOT.GIVEN, root.x = root, intermediates = FALSE, datc_se)
			return(ll)
		}
		attr(likSE, "cache") <- cache
		attr(likSE,"argn")=list(rates=cache$phy$edge[,2], root="root", SE="SE")
		class(likSE)=c("rbm","bm","function")
		return(likSE)
	} else { # unadjustable SE
		lik <- function(rates, root) {
			check.argn(rates, root)
			vv[mm] = rates
			ll = ll.bm.direct(pars = vv, root = ROOT.GIVEN, root.x = root, intermediates = FALSE, datc_init)
			return(ll)
		}
		attr(lik, "cache") <- cache
		attr(lik,"argn")=list(rates=cache$phy$edge[,2], root="root")
		class(lik)=c("rbm","bm","function")
		return(lik)
	}	
}

.cache.tree <- function (phy) 
{
	ordxx=function (children, is.tip, root) 
# from diversitree:::get.ordering
	{
		todo <- list(root)
		i <- root
		repeat {
			kids <- children[i, ]
			i <- kids[!is.tip[kids]]
			if (length(i) > 0) 
            todo <- c(todo, list(i))
			else break
		}
		as.vector(unlist(rev(todo)))
	}
	
    edge <- phy$edge
    edge.length <- phy$edge.length
    idx <- seq_len(max(edge))
    n.tip <- Ntip(phy)
    tips <- seq_len(n.tip)
    root <- n.tip + 1
    is.tip <- idx <= n.tip
	
	desc=.compile_descendants(phy)
    children <- desc$fdesc
	if(!max(sapply(children, length) == 2)){
		children=NULL
		order=NULL
		binary=FALSE
	} else {
		children <- rbind(matrix(NA, n.tip, 2), t(matrix(unlist(children), nrow=2)))
		order <- ordxx(children, is.tip, root)	
		binary=TRUE
	}
	
    len <- edge.length[mm<-match(idx, edge[, 2])]
	
	ans <- list(tip.label = phy$tip.label, node.label = phy$node.label, mm=mm,
				len = len, children = children, order = order, 
				root = root, n.tip = n.tip, n.node = phy$Nnode, tips = tips, 
				edge = edge, edge.length = edge.length, binary = binary, desc = desc)
    ans
}


.prepare.bm <- function(phy, dat, SE=NA) {

	## SE: can be array of mixed NA and numeric values -- where SE == NA, SE will be estimated (assuming a global parameter for all species) 
	# resolve estimation of SE
	if(is.null(SE)) SE=NA
	if(any(is.na(SE))) SE[is.na(SE)]=-666
	
	# check major issues
	if(all(is.null(names(dat)))) stop("'dat' must have names matchable to tip labels in 'phy'")
	if(!inherits(phy,"phylo")) stop("Supply 'phy' as a phylo object.")
	td=treedata(phy, dat, sort = TRUE, warnings=FALSE)
	if(any(sapply(td, length)==0)) stop("'dat' and 'phy' appear unmatchable.")
	if(ncol(td$data)>1) stop("Please supply one trait at a time")
	td$phy=reorder(td$phy,"pruningwise")
	
	# cache object
	cache=.cache.data.bm(td$phy, td$data[,1], SE)
	cache$nodes=cache$phy$edge[,2]
	return(cache)
}


.cache.data.bm <- function (phy, dat, SE=NULL) 
{
	cache=.cache.tree(phy)
	
	# RESOLVE SE
	# SE should be entirely numeric
    if (is.null(SE)) stop("'SE' must be a numeric value or vector")
	adjustSE=FALSE
	
	if(length(SE)==1) {
		SE=rep(SE, Ntip(phy))
		names(SE)=phy$tip.label
	}
	
	if(!is.null(names(SE))) {
		xx=setequal(names(SE), phy$tip.label)
		ss.match=match(names(SE), phy$tip.label)
		if(!xx) stop(paste("encountered error in constructing 'SE' with the following tips not found in the dataset:\n\t", paste(names(SE[is.na(ss.match)]), collapse="\n\t"), sep=""))
		SE=SE[ss.match]
	} else {
		stop("'SE' must be a named vector of measurement error")
	}
	
	tmp <- .check.states.quasse(phy, dat, SE)
    states <- tmp$states
    states.sd <- tmp$states.sd
    cache$ny <- 3
    y <- mapply(function(mean, sd) c(mean, sd * sd, 0), states, states.sd, SIMPLIFY = TRUE)

    cache$y <- .dt.tips.ordered(y, cache$tips, cache$len[cache$tips])
	cache$y$SE <- as.integer((SE==-666)[cache$y$target])
	
	cache$SE=SE
	cache$dat=dat
	cache$phy=phy
	
    cache
}

.proc.lnR=function(gen, subprop, cur.lnL, new.lnL, lnp, lnh, heat=1, control){
	if(is.infinite(cur.lnL)) stop("Starting point has exceptionally poor likelihood.")
	if(is.finite(new.lnL)) {
		lnLikelihoodRatio = new.lnL - cur.lnL
	} else {
		new.lnL=-Inf
		lnLikelihoodRatio = -Inf
	}

	if(control$sample.priors) lnLikelihoodRatio=0

	lnR=(heat * lnLikelihoodRatio) + (heat * lnp) + lnh

	r=.assess.lnR(lnR)	
	
	if(r$error) .error.rjmcmc(gen, subprop, cur.lnL, new.lnL, lnLikelihoodRatio, lnp, lnh, control$errorlog)
	return(r)
}


#general phylogenetic utility for determining whether to accept ('r'=1) or reject ('r'=0) a proposed model in Markov chain Monte Carlo sampling
#author: JM EASTMAN 2010
.assess.lnR <- function(lnR) {
	if(is.na(lnR) | is.infinite(abs(lnR))) {
		error=TRUE
		r=0
	} else {
		if(lnR < -20) {
			r=0 
		} else if(lnR >= 0) {
			r=1 
		} else {
			r=exp(lnR)
		}
		error=FALSE
	}
	return(list(r=r, error=error))
}



