## MAIN FUNCTION for OPTIMIZATION
# to replace geiger:::fitContinuous



fitContinuous=function(phy, dat, SE = NA, model=c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "white"), bounds=list(), control=list(method=c("SANN","L-BFGS-B"), niter=100, FAIL=1e200)){

	require(auteur)
	
	# SE: can be vector or single value (numeric, NULL, or NA); vector can include NA
	# opt: a list with elements 'method', 'niter', 'FAIL'; 'method' may include several optimization methods
	# bounds: a list with elements specifying constraint(s): e.g., bounds=list(alpha=c(0,1))
	# control: a list with elements specifying method to compute likelihood
	
	
	# CONTROL OBJECT for optimization
	ct=list(method=c("SANN","L-BFGS-B"), niter=ifelse(any(is.na(SE)), 100, 50), FAIL=1e200)
	if("method"%in%names(control)) control$method=match.arg(control$method, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), several.ok=TRUE)
	ct[names(control)]=control
	if(ct$niter<2) stop("'niter' must be equal to or greater than 2")
	
	model=match.arg(model, c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "white"))
	
	# CONTROL OBJECT for likelihood
	con=list(method="pruning",backend="C")
	con[names(control)]=control
	lik=bm.lik(phy,dat,SE,model)
	argn=unlist(argnames(lik))
	

	## CONSTRUCT BOUNDS ##
	mn=c(-500, -500, -3, -100, -500, -500, -500, -500)
	mx=c(100, 5, 3, 100, 0, 0, log(2.999999), 100)
	bnds=as.data.frame(cbind(mn, mx))
	bnds$typ=c("exp", "exp", "nat", "nat", "exp", "exp", "exp", "exp")
	rownames(bnds)=c("sigsq", "alpha", "a", "slope", "lambda", "kappa", "delta", "SE")
	bnds$model=c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "SE")
	
	# User bounds
	if(length(bounds)>0){
		mm=match(names(bounds), rownames(bnds))
		if(all(!is.na(mm))){
			for(i in 1:length(mm)){
				ww=mm[i]
				tmp=sort(bounds[[i]])
				if(bnds$typ[ww]=="exp") {
					if(any(tmp==0)) tmp[tmp==0]=exp(-500)
					bnd=log(tmp) 
				} else {
					bnd=tmp
				}
				bnds[ww,c("mn","mx")]=bnd
			}
		} else {
			stop(paste("'bounds' not expected:\n\t", paste(names(bounds)[is.na(mm)], sep="", collapse="\n\t"), sep=""))
		}
	}
		
	par=argn[1]

	
	## likelihood function for optimizer 
	xx=function(p){
		if(par=="sigsq") {
			if("SE"%in%argn){
				tmp=-lik(exp(p[1]), exp(p[2]))
			} else {
				tmp=-lik(exp(p[1]))
			}
		} else {
			if(bnds[par,"typ"]=="nat"){
				if("SE"%in%argn){
					tmp=-lik(p[1], exp(p[2]), exp(p[3]))
				} else {
					tmp=-lik(p[1], exp(p[2]))
				}
			} else {
				if("SE"%in%argn){
					tmp=-lik(exp(p[1]), exp(p[2]), exp(p[3]))
				} else {
					tmp=-lik(exp(p[1]), exp(p[2]))
				}
			}
		}
		if(is.infinite(tmp)) {
			tmp=ct$FAIL
		}
		tmp
	}
	
	# boxconstrain from diversitree
	boxconstrain=function (f, lower, upper, fail.value = FAIL) 
	{
		function(x) {
			if (any(x < lower | x > upper)) fail.value else f(x)
		}
	}
	f=boxconstrain(xx, min, max, fail.value=ct$FAIL)

	
	## STARTING POINT -- problematic models ##
	if(par%in%c("alpha","lambda","delta","kappa")){
		bmstart=try(.bm.smartstart(phy,dat),silent=TRUE)
		if(inherits(bmstart, "try-error")) bmstart=0.01 
		if(par=="alpha") oustart=.ou.smartstart(bmstart,var(dat))
	}

	## OPTIMIZATION ##
	mm=matrix(NA, nrow=ct$niter, ncol=length(argn)+2)
	mt=character(ct$niter)
	
	# 'method' optimization
	for(i in 1:ct$niter){
		bnds$st=sapply(1:nrow(bnds), function(x) runif(1, bnds$mn[x], bnds$mx[x]))
		start=bnds[argn,"st"]
		
		## OU ##
		if(par=="alpha"){
			if(i==1 | runif(1)<0.25) start[match(c("sigsq", par), argn)]=c(bmstart, oustart)
			if(runif(1) < 0.5) start[match(c("sigsq", par), argn)]=c(0,oustart)
		}
		
		## PAGEL MODELS ##
		if(par%in%c("lambda","delta","kappa")){
			ww=match(par, rownames(bnds))
			if(runif(1)<0.5){
				if(runif(1)<0.5){
					start[match(c("sigsq", par), argn)]=c(bmstart, bnds$mx[ww])
				} else {
					start[match(c("sigsq", par), argn)]=c(bmstart, bnds$mn[ww])
				}
			}
		}
		
		## WHITE NOISE ##
		if(par=="white"){
			if(runif(1)<0.5){
				start[match("sigsq", argn)]=var(dat)
			}
		}
				
		names(start)=argn
		min=bnds[argn,"mn"]
		max=bnds[argn,"mx"]
		typs=bnds[argn, "typ"]

		
		# resolve method	
		if(length(argn)==1) {
			method="Brent" 
		} else {
			if(runif(1)<ifelse(model=="OU", 0.25, 0.5)){
				method="subplex"
			} else {
				method=sample(ct$method,1)
			}
		}
		
		if(method=="subplex"){
			op=try(suppressWarnings(subplex(par=start, fn=f, control=list(reltol = .Machine$double.eps^0.25, parscale = rep(0.1, length(argn))))),silent=TRUE)
		} else {
			op=try(suppressWarnings(optim(par=start, fn=f, upper=max, lower=min, method=method)),silent=TRUE)
		}
		op$method=method

		if(!inherits(op,"try-error")){
			op$value=-op$value
			names(op)[names(op)=="value"]="lnL"
			names(op$par)=argn
			
			op$par=sapply(1:length(typs), function(x) if(typs[x]=="exp") return(exp(op$par[x])) else return(op$par[x]))
			mm[i, ]=c(op$par, op$lnL, op$convergence)
			mt[i]=op$method
		}
	}
	res=mm
	colnames(res)=c(argn,"lnL","convergence")
	rownames(res)=mt

	## HANDLE OPTIMIZER OUTPUT ##
	colnames(mm)=c(argn, "lnL", "convergence")
	conv=mm[,"convergence"]==0
	mm=mm[,-which(colnames(mm)=="convergence")]
	valid=apply(mm, 1, function(x) !any(is.na(x)))
	if(sum(valid & conv)>=1){
		mm=matrix(mm[valid,], nrow=sum(valid), dimnames=dimnames(mm))
		mt=mt[valid]
		mm=mm[z<-min(which(mm[,"lnL"]==max(mm[,"lnL"]))),]
	} else {
		z=NA
		mm=c(rep(NA, length(argn)), -Inf)
		names(mm)=c(argn,"lnL")
	}
	mm=as.list(mm)
	mm$method=ifelse(is.na(z), NA, mt[z])
	mm$k=length(argn)+1
	
	# check estimates against bounds #
	range=as.data.frame(cbind(min, max))
	range$typ=typs
	range$mn=ifelse(range$typ=="exp", exp(range$min), range$min)
	range$mx=ifelse(range$typ=="exp", exp(range$max), range$max)
	par=mm[argn]
	rownames(range)=argn
	chk=sapply(1:length(par), function(idx){
		   p=par[[idx]]
		   if(!is.na(p)){
				return((p<=range$mn[idx] | p>=range$mx[idx]))
		   } else {
				return(FALSE)
		   }
	})
	if(any(chk)){
			warning(paste("Parameter estimates appear at bounds:\n\t", paste(names(par)[chk], collapse="\n\t", sep=""), sep=""))
	}
	
	mm=.aic(mm, n=length(dat))
	return(list(lik=lik, bnd=range[,c("mn", "mx")], res=res, opt=mm))
}


## WORKHORSE -- built from diversitree:::make.bm by tricking models into multivariate normal 
# likelihood function creation
bm.lik<-function (phy, dat, SE = NA, model=c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "white")) 
{
	## SE: can be array of mixed NA and numeric values -- where SE == NA, SE will be estimated (assuming a global parameter for all species) 

	model=match.arg(model, c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "white"))

	# resolve binary tree
	if(!is.binary.tree(phy)){
		phy=multi2di(phy)
	}
	phy=reorder(phy)
	
	# resolve estimation of SE
	if(is.null(SE)) SE=NA
	if(any(is.na(SE))) {
		if(model=="white"){
			warning("'SE' and 'sigsq' are confounded in the 'white' model:\n\t'SE' will be ignored")
			adjSE=FALSE 
			SE[is.na(SE)]=0
		} else {
			adjSE=TRUE 
			SE[is.na(SE)]=0
		}
	} else {
		adjSE=FALSE
	}
	
	# cache object for likelihood computation
	cache = diversitree:::make.cache.bm(phy, dat, SE, control=list(method="pruning", backend="C"))
	cache$ordering=attributes(cache$info$phy)$order
		
	FUN=switch(model, 
			   BM=.null.cache(cache),
			   OU=.ou.cache(cache),
			   EB=.eb.cache(cache),
			   trend=.trend.cache(cache),
			   lambda=.lambda.cache(cache),
			   kappa=.kappa.cache(cache),
			   delta=.delta.cache(cache),
			   white=.white.cache(cache)
			   )

	cache$adjustSE=ifelse(cache$y$y[2,]==0, TRUE, FALSE)
	
    z = length(cache$len)
    rr = numeric(z)
    rootidx = as.integer(cache$root)
	nn = length(cache$len)
	
    datc = list(intorder = as.integer(cache$order[-length(cache$order)]), 
				tiporder = as.integer(cache$y$target), 
				root = rootidx, 
				y = as.numeric(cache$y$y[1, ]), 
				n = as.integer(z), 
				descRight = as.integer(cache$children[ ,1]), 
				descLeft = as.integer(cache$children[, 2]))
	
    ll.bm.direct = function(q, sigsq, se, root = ROOT.MAX, root.x = NA, intermediates = FALSE, datc) {
		
		if(is.null(argnames(FUN))) new=FUN() else new=FUN(q)
		datc$len=as.numeric(new$len)
		.xxSE=function(cache){
			vv=cache$y$y[2,]
			ff=function(x){
				if(adjSE){
					vv[which(cache$adjustSE)]=x^2 
					return(vv)
				} else {
					return(vv)
				}
			}
			return(ff)
		}
		modSE=.xxSE(cache)
		
		datc$var=as.numeric(modSE(se))
		
        out = .Call("bm_direct", dat = datc, pars = as.numeric(rep(sigsq, nn)), package = "auteur")
        vals = c(out$initM[rootidx], out$initV[rootidx], out$lq[rootidx])
        loglik <- auteur:::.root.bm.direct(vals, out$lq[-rootidx], root, root.x)
        if (intermediates) {
            attr(loglik, "intermediates") <- intermediates
            attr(loglik, "vals") <- vals
        }
		if(is.na(loglik)) loglik=-Inf
        return(loglik)
    }
    class(ll.bm.direct) <- c("bm", "dtlik", "function")
	
	
	## EXPORT LIKELIHOOD FUNCTION
	if(is.null(argnames(FUN))){
		
		if(adjSE){
			lik <- function(sigsq, SE) {
				ll = ll.bm.direct(q=NULL, sigsq = sigsq, se=SE, root = ROOT.MAX, root.x = NA, intermediates = FALSE, datc)
				return(ll)
			}
			attr(lik, "argnames") = list(sigsq="sigsq", SE="SE")
		} else {
			lik <- function(sigsq) {
				ll = ll.bm.direct(q=NULL, sigsq = sigsq, se=NULL, root = ROOT.MAX, root.x = NA, intermediates = FALSE, datc)
				return(ll)
			}
			attr(lik, "argnames") = list(sigsq="sigsq")
		}
		
	} else {
		if(adjSE){
			lik <- function(q, sigsq, SE) {
				ll = ll.bm.direct(q, sigsq = sigsq, se=SE, root = ROOT.MAX, root.x = NA, intermediates = FALSE, datc)
				return(ll)
			}
			attr(lik, "argnames") = list(par=argnames(FUN), sigsq="sigsq", SE="SE")			
		} else {
			lik <- function(q, sigsq) {
				ll = ll.bm.direct(q, sigsq = sigsq, se=NULL, root = ROOT.MAX, root.x = NA, intermediates = FALSE, datc)
				return(ll)
			}
			attr(lik, "argnames") = list(par=argnames(FUN), sigsq="sigsq")			
		}
	}
	
    attr(lik, "cache") <- cache
    class(lik) = c("bm", "function")
    lik
}



