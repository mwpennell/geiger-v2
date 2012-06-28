## MAIN FUNCTION for OPTIMIZATION
# to replace geiger:::fitContinuous

fitContinuous=function(
	phy, 
	dat, 
	SE = NA, 
	model=c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "white"), 
	bounds=list(), 
	control=list(method=c("SANN","L-BFGS-B"), niter=100, FAIL=1e200, hessian=FALSE, hessian_P=0.05)
	)
{

#	require(auteur)
	
	# SE: can be vector or single value (numeric, NULL, or NA); vector can include NA
	# opt: a list with elements 'method', 'niter', 'FAIL'; 'method' may include several optimization methods
	# bounds: a list with elements specifying constraint(s): e.g., bounds=list(alpha=c(0,1))
	# control: a list with elements specifying method to compute likelihood
	
	# data matching
	td=treedata(phy, dat)
	phy=td$phy
	dat=td$data
	dd=dim(dat)
	trts=dd[2]
	if(trts>1){
		nm=colnames(dat)
		res=lapply(1:trts, function(idx){
			fitContinuous(phy, dat[,idx], SE=SE, model=model, bounds=bounds, control=control)	   
		})
		names(res)=nm
		return(res)
	} else {
		dat=dat[,1]
	}

	
	# CONTROL OBJECT for optimization
	ct=list(method=c("SANN","L-BFGS-B"), niter=ifelse(any(is.na(SE)), 100, 50), FAIL=1e200, hessian=FALSE)
	if("method"%in%names(control)) control$method=match.arg(control$method, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), several.ok=TRUE)
	ct[names(control)]=control
	if(ct$niter<2) stop("'niter' must be equal to or greater than 2")
	
	model=match.arg(model, c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "white"))
	
	# CONTROL OBJECT for likelihood
	con=list(method="pruning",backend="C")
	con[names(control)]=control
	lik=bm.lik(phy,dat,SE,model)
	argn=argnames(lik)
	

	## CONSTRUCT BOUNDS ##
	mn=c(-500, -500, -5, -100, -500, -500, -500, -500)
	mx=c(100, 5, 0.1, 100, 0, 0, log(2.999999), 100)
	bnds=as.data.frame(cbind(mn, mx))
	bnds$typ=c("exp", "exp", "nat", "nat", "exp", "exp", "exp", "exp")
	rownames(bnds)=c("sigsq", "alpha", "a", "slope", "lambda", "kappa", "delta", "SE")
	bnds$model=c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "SE")
	typs=bnds[argn, "typ"]
	
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

	
	## likelihood function for optimizer (with modified space)
	xx=function(p){
		pars=ifelse(typs=="exp", exp(p), p)
		tmp=-lik(pars)
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
	out=list()
	
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
			op<-out[[i]]<-try(suppressWarnings(subplex(par=start, fn=f, control=list(reltol = .Machine$double.eps^0.25, parscale = rep(0.1, length(argn))), hessian=ct$hessian)),silent=TRUE)
		} else {
			op<-out[[i]]<-try(suppressWarnings(optim(par=start, fn=f, upper=max, lower=min, method=method, hessian=ct$hessian)),silent=TRUE)
		}
		op$method=method

		if(!inherits(op,"try-error")){
			op$value=-op$value
			names(op)[names(op)=="value"]="lnL"
			names(op$par)=argn
			
			op$par=sapply(1:length(typs), function(x) if(typs[x]=="exp") return(exp(op$par[x])) else return(op$par[x]))
			mm[i, ]=c(op$par, op$lnL, op$convergence)
			mt[i]=op$method
		} else {
			mt[i]="FAIL"
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
		out=out[valid]
		mm=mm[z<-min(which(mm[,"lnL"]==max(mm[,"lnL"]))),]
	} else {
		z=NA
		mm=c(rep(NA, length(argn)), -Inf)
		names(mm)=c(argn,"lnL")
	}
	zz=mm[-which(names(mm)%in%c("lnL"))]
	mm=as.list(mm)
	mm$method=ifelse(is.na(z), NA, mt[z])
	mm$k=length(argn)+1
	
	## HESSIAN-based CI of par estimates
	if(ct$hessian){
		hessian=out[[z]]$hessian
		CI=.bnd.hessian(hessian, zz, typs, ct$hessian_P)
		if(!all(is.na(CI))){
			if(is.constrained(lik)){
				CI=rbind(lik(CI[1,], pars.only=TRUE), rbind(lik(CI[2,], pars.only=TRUE)))
			} 
			dimnames(hessian)=NULL
			rownames(CI)=c("lb", "ub")
		}
	} else {
		hessian=NULL
		CI=NULL
	}
	
	
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

	# RETURN OBJECT
	mm$CI=CI
	mm$hessian=hessian
#	if(is.null(CI)){
		return(list(lik=lik, bnd=range[,c("mn", "mx")], res=res, opt=mm))	
		
#	} else {
#		return(list(lik=lik, bnd=range[,c("mn", "mx")], res=res, opt=mm, CI=CI))	
#		
#	}
	
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
			warning("'SE' and 'sigsq' are confounded in the 'white' model:\n\t'SE' will be set to nil")
			adjSE=FALSE 
			SE[is.na(SE)]=0
		} else {
			adjSE=TRUE 
			SE[is.na(SE)]=-666
		}
	} else {
		adjSE=FALSE
	}
	
	# cache object for likelihood computation
	cache = .make.cache.bm(phy, dat, SE, control=list(method="pruning", backend="C"))
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

	cache$adjustSE=ifelse(cache$y$y[2,]==-666, TRUE, FALSE)
	
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
        loglik <- .root.bm.direct(vals, out$lq[-rootidx], root, root.x)
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
			attb=c("sigsq", "SE")
			lik <- function(pars) {
				pars=.repars(pars, attb)
				ll = ll.bm.direct(q=NULL, sigsq = pars[1], se=pars[2], root = ROOT.MAX, root.x = NA, intermediates = FALSE, datc)
				return(ll)
			}
			attr(lik, "argnames") = attb
		} else {
			attb="sigsq"
			lik <- function(pars) {
				pars=.repars(pars, attb)
				ll = ll.bm.direct(q=NULL, sigsq = pars[1], se=NULL, root = ROOT.MAX, root.x = NA, intermediates = FALSE, datc)
				return(ll)
			}
			attr(lik, "argnames") = attb
		}
		
	} else {
		if(adjSE){
			attb=c(argnames(FUN), "sigsq", "SE")
			lik <- function(pars) {
				pars=.repars(pars, attb)
				ll = ll.bm.direct(q=pars[1], sigsq = pars[2], se=pars[3], root = ROOT.MAX, root.x = NA, intermediates = FALSE, datc)
				return(ll)
			}
			attr(lik, "argnames") = attb			
		} else {
			attb=c(argnames(FUN), "sigsq")
			lik <- function(pars) {
				pars=.repars(pars, attb)
				ll = ll.bm.direct(pars[1], sigsq = pars[2], se=NULL, root = ROOT.MAX, root.x = NA, intermediates = FALSE, datc)
				return(ll)
			}
			attr(lik, "argnames") = attb			
		}
	}
	
    attr(lik, "cache") <- cache
    class(lik) = c("bm", "function")
    lik
}



