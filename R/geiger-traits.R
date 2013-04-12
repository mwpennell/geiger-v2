## MAIN FUNCTION for OPTIMIZATION
# to replace geiger:::fitContinuous

fitContinuous=function(
phy, 
dat, 
SE = 0,
model=c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "drift", "white"), 
bounds=list(), 
control=list(method=c("SANN","L-BFGS-B"), niter=100, FAIL=1e200, hessian=FALSE, CI=0.95), ...
)
{
		
# SE: can be vector or single value (numeric, NULL, or NA); vector can include NA
# opt: a list with elements 'method', 'niter', 'FAIL'; 'method' may include several optimization methods
# bounds: a list with elements specifying constraint(s): e.g., bounds=list(alpha=c(0,1))
# control: a list with elements specifying method to compute likelihood
# ...: node state dataframe
	
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
        class(res)=c("gfits", class(res))
		return(res)
	} else {
		dat=dat[,1]
	}
	
	
# CONTROL OBJECT for optimization
	ct=list(method=c("SANN","L-BFGS-B"), niter=100, FAIL=1e200, hessian=FALSE, CI=0.95)
	if(any(!names(control)%in%names(ct)->tmp)) warning("Unexpected 'control' parameters:\n\t", paste(names(control)[tmp], collapse="\n\t"), sep="")
	control=control[which(!tmp)]
	if("method"%in%names(control)) control$method=match.arg(control$method, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), several.ok=TRUE)
	ct[names(control)]=control
	if(ct$niter<2) stop("'niter' must be equal to or greater than 2")
	ct$hessian_P=1-ct$CI
	
	model=match.arg(model, c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "drift", "white"))
	
# CONTROL OBJECT for likelihood
	con=list(method="pruning",backend="C")
	con[names(control)]=control
	lik=bm.lik(phy,dat,SE,model,...)
	argn=argn(lik)
	
	
## CONSTRUCT BOUNDS ##
	mn=c(-500, -500, -1, -100, -100, -500, -500, -500, -500, -1e6)
	mx=c(100, 1, 1, 100, 100, 0, 0, log(2.999999), 100, 1e6)
	bnds=as.data.frame(cbind(mn, mx))
	bnds$typ=c("exp", "exp", "nat", "nat", "nat", "exp", "exp", "exp", "exp", "nat")
	rownames(bnds)=c("sigsq", "alpha", "a", "drift", "slope", "lambda", "kappa", "delta", "SE", "z0")
	bnds$model=c("BM", "OU", "EB", "drift", "trend", "lambda", "kappa", "delta", "SE", "z0")
	typs=bnds[argn, "typ"]
	
# User bounds
	if(length(bounds)>0){
		mm=match(names(bounds), rownames(bnds))
		if(any(is.na(mm))){
			warning("Unexpected 'bounds' parameters:\n\t", paste(names(bounds)[is.na(mm)], collapse="\n\t"), sep="")
		}
		mm=mm[!is.na(mm)]
		
		if(length(mm)){
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
		} 
	}
	
	par=argn[1]
	
	
## likelihood function for optimizer (with modified space)
	xx=function(p){
		pars=ifelse(typs=="exp", exp(p), p)
		tmp=-lik(pars)
		if(is.infinite(tmp)) tmp=ct$FAIL
        if(is.na(tmp)) tmp=ct$FAIL
		tmp
	}
	
# boxconstrain (from diversitree) to avoid out-of-bounds values
	boxconstrain=function (f, lower, upper, fail.value) 
	{
		function(x) {
			if (any(x < lower | x > upper)) fail.value else f(x)
		}
	}
	f=boxconstrain(xx, bnds[argn,"mn"], bnds[argn,"mx"], fail.value=ct$FAIL)
	
	
## STARTING POINT -- problematic models ##
	if(par%in%c("alpha","lambda","delta","kappa")){
		bmstart=try(.bm.smartstart(phy,dat),silent=TRUE)
		if(inherits(bmstart, "try-error")) bmstart=0.01
        bmstart=log(bmstart) 
	}
    rt=max(abs(dat))
	
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
            oustart=log(.ou.smartstart(dat, unlist(exp(bnds["alpha",c("mx","mn")]))))
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
        start[["z0"]]=runif(1, -rt, rt)
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
		
		if(!inherits(op,"try-error")){
            op$method=method
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
	mm$k=length(argn)
	
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
    res=list(lik=lik, bnd=range[,c("mn", "mx")], res=res, opt=mm)
    class(res)=c("gfit", class(res))
    return(res)
}


## create 'standard' constraint likelihood function given a patterned matrix
constrain.k=function(f, model=c("ER", "SYM", "ARD", "meristic"), ...){
	e=environment(f)
	k=e$k
	
	m=matrix(1:k^2,nrow=k,byrow=TRUE)
	model=match.arg(model,c("ER", "SYM", "ARD", "meristic"))
	mm=switch(model,
			  ER=.er.matrix(m),
			  SYM=.sym.matrix(m),
			  ARD=.ard.matrix(m),
			  meristic=.meristic.matrix(m, ...)
			  )
	attributes(mm)
	ff=constrain.m(f,mm)
	attr(ff,"constraint.m")=.reshape.constraint.m(mm)
	ff
}

.reshape.constraint.m=function(m){
	k=unique(dim(m))
	if(length(k)>1) stop("'m' must be a square matrix")
	diag(m)=NA
	tt=table(m)
	map=cbind(as.integer(names(tt)), seq_along(1:max(as.integer(length(tt)))))
	z=match(m, map[,1])
	m[]=map[z,2]
	class(m)=c("constraint.m",class(m))
	m
}

print.constraint.m=function(x, printlen=3, ...){
	cat("matrix representation of unique and shared transitions \n")
	tt=table(x)
	if(any(tt==1)){
		cat("\tunique transition classes:", paste(sort(as.integer(names(tt[tt==1]))), collapse=", "), "\n")
	}
	if(any(tt>1)){
		cat("\tshared transition classes:", paste(sort(as.integer(names(tt[tt>1]))), collapse=", "),"\n")
	}
	cat("\n")
	class(x)="matrix"
	print(x)
}


## MAIN FUNCTION for OPTIMIZATION
# to replace geiger:::fitContinuous

fitDiscrete=function(	
phy, 
dat, 
model=c("ER","SYM","ARD","meristic"),
transform=c("none", "EB","lambda", "kappa", "delta", "white"), 
bounds=list(), 
control=list(method=c("SANN","L-BFGS-B"), niter=100, FAIL=1e200, hessian=FALSE, CI=0.95),
...)
{
	
## NOTE: 'model' can be a constraint matrix
#	
#		transform="none"
#		model="SYM" 
#		bounds=list()
#		control=list(hessian=TRUE) 
	
	td=treedata(phy, dat)
	phy=td$phy
	dat=td$data
	dd=dim(dat)
	trts=dd[2]
	if(trts>1){
		nm=colnames(dat)
		res=lapply(1:trts, function(idx){
				   fitDiscrete(phy, dat[,idx], model=model, transform=transform, bounds=bounds, control=control, ...)	   
				   })
		names(res)=nm
        class(res)=c("gfits", class(res))
		return(res)
	} else {
		tmp=dat[,1]
		names(tmp)=rownames(dat)
		dat=tmp
		#if(!all(is.integer(dat))) stop("supply 'dat' as a vector (or matrix) of positive integers")
		
	}
	
	constrain=model
	model=match.arg(transform, c("none", "EB", "lambda", "kappa", "delta", "white"))
	
# CONTROL OBJECT
	ct=list(method=c("SANN","L-BFGS-B"), niter=ifelse(model=="none", 50, 100), FAIL=1e200, hessian=FALSE, CI=0.95)
	if(any(!names(control)%in%names(ct)->tmp)) warning("Unexpected 'control' parameters:\n\t", paste(names(control)[tmp], collapse="\n\t"), sep="")
	control=control[which(!tmp)]
	if("method"%in%names(control)) control$method=match.arg(control$method, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), several.ok=TRUE)
	ct[names(control)]=control
	if(ct$niter<2) stop("'niter' must be equal to or greater than 2")
	ct$hessian_P=1-ct$CI
	
	
	lik=mkn.lik(phy, dat, constrain=constrain, transform=model, control=list(method="exp", root=ROOT.OBS), ...)
	if(model=="white") return(list(opt=lik))
	argn=unlist(argn(lik))
	
## CONSTRUCT BOUNDS ##
	mn=c(-10, -500, -500, -5, -500)
	mx=c(10, 0, 0, log(2.999999), log(1000))
	bnds=as.data.frame(cbind(mn, mx))
	bnds$typ=c("nat", "exp", "exp", "exp", "exp")
	rownames(bnds)=c("a", "lambda", "kappa", "delta", "trns")
	bnds$model=c("EB", "lambda", "kappa", "delta", "trns")
	
	parnm=ifelse(argn%in%rownames(bnds), argn, "trns")
	typs=bnds[parnm, "typ"]
	
# User bounds
	if(length(bounds)>0){
		mm=match(names(bounds), rownames(bnds))
		if(any(is.na(mm))){
			warning("Unexpected 'bounds' parameters:\n\t", paste(names(bounds)[is.na(mm)], collapse="\n\t"), sep="")
		}
		mm=mm[!is.na(mm)]
		if(length(mm)){
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
		} 
	}
	
	par=argn[1]
	
## likelihood function for optimizer (with modified space)
	xx=function(p){
		pars=ifelse(typs=="exp", exp(p), p)
		tmp=-lik(pars)
		if(is.infinite(tmp)) tmp=ct$FAIL
        if(is.na(tmp)) tmp=ct$FAIL
		tmp
	}
	
	
	
# boxconstrain from diversitree
	boxconstrain=function (f, lower, upper, fail.value) 
	{
		function(x) {
			if (any(x < lower | x > upper)) fail.value else f(x)
		}
	}
	f=boxconstrain(xx, min, max, fail.value=ct$FAIL)
	
# transition rates
	smarttrns=log(seq(max(c(exp(bnds["trns","mn"]),0.01)), min(c(exp(bnds["trns","mx"]),2)), length.out=10))
	
## OPTIMIZATION ##
	mm=matrix(NA, nrow=ct$niter, ncol=length(argn)+2)
	mt=character(ct$niter)
	out=list()
	
# 'method' optimization
	for(i in 1:ct$niter){
		bnds$st=sapply(1:nrow(bnds), function(x) runif(1, bnds$mn[x], bnds$mx[x]))
		start=bnds[parnm,"st"]
		
## PAGEL MODELS ##
		if(par%in%c("lambda","delta","kappa")){
			ww=match(par, rownames(bnds))
			if(runif(1)<0.5){
				if(runif(1)<0.5){
					start[match(par, argn)]=bnds$mx[ww]
				} else {
					start[match(par, argn)]=bnds$mn[ww]
				}
			}
		}
		
## TRANSITION RATES
		if(runif(1)<0.5){
			start[which(parnm=="trns")][]=smarttrns[sample(1:length(smarttrns),1)]
		}
		
		names(start)=argn
		min=bnds[parnm,"mn"]
		max=bnds[parnm,"mx"]
		
# resolve method	
		if(length(argn)==1) {
			method="Brent" 
		} else {
#			if(runif(1)<0.5){  ### FIXME ###
#				method="subplex"
#			} else {
				method=sample(ct$method,1)
#			}
		}
		
		if(method=="subplex"){
			op<-out[[i]]<-try(suppressWarnings(subplex(par=start, fn=f, control=list(reltol = .Machine$double.eps^.25, parscale = rep(0.1, length(argn))), hessian=ct$hessian)),silent=TRUE)
		} else {
			op<-out[[i]]<-try(suppressWarnings(optim(par=start, fn=f, upper=max, lower=min, method=method, hessian=ct$hessian)),silent=TRUE)
		}
		
		if(!inherits(op,"try-error")){
            op$method=method
			op$value=-op$value
			names(op)[names(op)=="value"]="lnL"
			names(op$par)=argn
			
			op$par=sapply(1:length(typs), function(x) if(typs[x]=="exp") return(exp(op$par[x])) else return(op$par[x]))
			mm[i, ]=c(op$par, op$lnL, op$convergence)
			mt[i]=op$method
		} else {
            print(start)
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
		mod=mt[z]
	} else {
		mod=NA
		z=NA
		mm=c(rep(NA, length(argn)), -Inf)
		names(mm)=c(argn,"lnL")
	}
	k=length(argn)
	llx=which(names(mm)=="lnL")
	ll=mm[[llx]]
	mm=mm[-llx]
	
## HESSIAN-based CI of par estimates
	if(ct$hessian){
#		print(mm)
#		print(partp)
#		print(out[[z]]$hessian)
		hessian=out[[z]]$hessian
		CI=.bnd.hessian(hessian, mm, typs, ct$hessian_P)
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
	
# resolve all transition estimates (if constrained model)
	trn=!names(mm)%in%c("a", "lambda", "kappa", "delta")
	if(is.constrained(lik)){
		constr=TRUE
		allq=lik(mm, pars.only=TRUE) 
	} else {
		constr=FALSE
		allq=mm[argn(lik)]
	}
	allq=allq[!names(allq)%in%names(mm)[!trn]]
	
	qparnm=rep("trns", length(allq))
	allparnm=c(parnm[!trn],qparnm)
	min=bnds[allparnm,"mn"]
	max=bnds[allparnm,"mx"]
	typs=bnds[allparnm, "typ"]
	argn=c(argn[!trn], names(allq))
	mmcomplete=c(mm[!trn], allq, ll)
	names(mmcomplete)=c(argn, "lnL")
	
	mm=as.list(mmcomplete)
	mm$method=mod
	mm$k=k
	mm=.aic(mm, n=length(dat))
	
	
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
	
# RETURN OBJECT
	mm$CI=CI
	mm$hessian=hessian
	
	res=list(lik=lik, bnd=range[,c("mn", "mx")], res=res, opt=mm)
    class(res)=c("gfit", class(res))
    return(res)
}

## WORKHORSE -- built from diversitree:::make.mkn  
# likelihood function creation
mkn.lik=function(
phy, 
dat, 
constrain=c("ER","SYM","ARD","meristic"),
transform=c("none", "EB", "lambda", "kappa", "delta", "white"), 
control=list(root=ROOT.OBS),
...)
{
    phy=reorder(phy, "postorder")
	
# control object for make.mkn()
	ct = list(method="exp", root=ROOT.OBS)
	ct[names(control)]=control
	if(ct$method!="exp") stop(paste("method",sQuote(ct$method),"is not currently supported",sep=" "))
	
# primary cache
	k<-nlevels(as.factor(dat))
    if(is.character(dat)) dat=structure(as.factor(dat), names=names(dat))
    if(is.factor(dat)){
        levels=levels(dat)
        dat=structure(as.integer(dat), names=names(dat))
    } else {
        levels=NULL
    }
    if(k==2) if(all(constrain=="SYM")) constrain="ER"
    control <- .check.control.mkn(ct, k)
    cache <- .make.cache.mkn(phy, dat, k, strict=TRUE, control=ct)
	cache$ordering=attributes(cache$info$phy)$order
	
# tree transforms
	trns=match.arg(transform, c("none", "EB","lambda", "kappa", "delta", "white"))
	
	FUN=switch(trns, 
			   none=.null.cache(cache),
			   EB=.eb.cache(cache),
			   lambda=.lambda.cache(cache),
			   kappa=.kappa.cache(cache),
			   delta=.delta.cache(cache), 
			   white=white.mkn(cache$states))
	
	if(trns=="white") return(FUN)	
	
	ll.mkn=function(cache, control) {
		k <- cache$info$k
		f.pars <- .make.pars.mkn(k)
		f.pij <- .make.pij.mkn(cache$info, control)
		idx.tip <- cache$idx.tip
		n.tip <- cache$n.tip
		n <- length(cache$len)
		map <- t(sapply(1:k, function(i) (1:k) + (i - 1) * k))
		idx.tip <- cbind(c(map[cache$states, ]), rep(seq_len(n.tip), k))
		children.C <- .toC.int(t(cache$children))
		order.C <- .toC.int(cache$order)
		
		.ll.mkn.exp=function(q, pars, intermediates=FALSE, preset = NULL) { # based on diversitree:::make.all.branches.mkn.exp
			if(is.null(argn(FUN))) new=FUN() else new=FUN(q)
			
			len.uniq <- sort(unique(new$len))
			len.idx <- match(new$len, len.uniq)
			
			if (!is.null(preset)) stop("Preset values not allowed")
			pij <- f.pij(len.uniq, pars)[, len.idx]
			lq <- numeric(n)
			branch.init <- branch.base <- matrix(NA, k, n)
			storage.mode(branch.init) <- "numeric"
			ans <- matrix(pij[idx.tip], n.tip, k)
			q <- rowSums(ans)
			branch.base[, seq_len(n.tip)] <- t.default(ans/q)
			lq[seq_len(n.tip)] <- log(q)
			ans <- .C("r_mkn_core", k = as.integer(k), n = length(order.C) - 
					  1L, order = order.C, children = children.C, pij = pij, 
					  init = branch.init, base = branch.base, lq = lq, 
					  NAOK = TRUE, DUP = FALSE, PACKAGE="geiger")
			
			
			list(init = ans$init, base = ans$base, lq = ans$lq, vals = ans$init[, cache$root], pij = pij)
		}
		
    # build likelihood function
		attb=c(argn(FUN), cache$info$argn)
		if(is.null(argn(FUN))){ # NO TRANSFORM
			ll=function(pars){
				qmat=f.pars(pars)
				ans=.ll.mkn.exp(q=NULL, pars=qmat, intermediates=FALSE)
				.rootfunc.mkn(ans, qmat, ct$root, NULL, intermediates=FALSE)
			}
		} else {
			ll=function(pars){ # TREE TRANSFORM
				qmat=f.pars(pars[-1])
				ans=.ll.mkn.exp(q=pars[1], pars=qmat, intermediates=FALSE)
				.rootfunc.mkn(ans, qmat, ct$root, NULL, intermediates=FALSE)
			}
			
		}
		class(ll) <- c("mkn", "dtlik", "function")
		attr(ll,"argn") <- attb
        if(!is.null(levels)) attr(ll, "levels")=levels
		return(ll)
	}
	
	tmp=ll.mkn(cache, control)
	
    ## CONSTRAINTS
	if(!all(constrain=="ARD")){
		if(is.character(constrain)){
			cc=match.arg(constrain, c("ER","SYM","ARD","meristic"))
			tmp=constrain.k(tmp, model=cc, ...)
		} else {
			if(is.matrix(constrain)){
				if(ncol(constrain)==max(dat)){
					tmp=constrain.m(tmp, m=constrain)
				}
			} else {
				stop("'constrain' must be supplied as a dummy matrix representing constraints on transition classes")
			}
		}
	}
	lik=function(pars, ...){
		pars=.repars(pars, argn(tmp))
		tmp(pars, ...)
	}
	attributes(lik)<-attributes(tmp)
	lik
}

aov.phylo=function(formula, phy, nsim=1000, test=c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"), ...){
    xx=lapply(all.vars(formula), get)
    flag="'formula' must be of the form 'dat~group', where 'group' is a named factor vector and 'dat' is a data matrix or named vector"
    
    if(!is.factor(xx[[2]])) stop(flag)
    if(is.null(names(xx[[2]]))) stop(flag)
    
    yy=merge(xx[[1]], xx[[2]], by=0)
    if(nrow(yy)==0) stop(flag)
    rownames(yy)=yy[,1]
    yy=yy[,-1]
    
    tmp<-treedata(phy, yy, sort=TRUE)
    phy=tmp$phy
    yy=yy[phy$tip.label,]
    
    group=structure(yy[,ncol(yy)], names=rownames(yy))
    dat=as.matrix(yy[,-ncol(yy)])
    rownames(dat)=rownames(yy)
        
    s<-vcv(phy, dat)
  
    multivar=ifelse(ncol(dat)>1, TRUE, FALSE)
    if(multivar){
        test=match.arg(test, c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"))
        m=summary.manova(mod<-manova(dat~group), test=test)
        f.data=m[[4]][1,2]
        FUN=function(xx) summary.manova(manova(as.matrix(xx)~group), test=test)[[4]][1,2]
        sims<-sim.char(phy, s, nsim=nsim)
        f.null<-apply(sims, 3, FUN)
        out=as.data.frame(m[[4]])
        attr(out, "heading")=c("Multivariate Analysis of Variance Table\n","Response: dat")
    } else {
        test=NULL
        m=anova(mod<-lm(dat~group))
        f.data<-m[1,4]
        FUN=function(xx) anova(lm(xx~group))[1,4]
        out=as.data.frame(m)
    }
    colnames(out)=gsub(" ", "-", colnames(out))
    sims<-sim.char(phy, s, nsim=nsim)
    f.null<-apply(sims, 3, FUN)
    p.phylo=(sum(f.null>f.data)+1)/(nsim+1)

    out$'Pr(phy)'=c(p.phylo, NA)
    print.anova(out, ...)
    attr(mod, "summary")=out
    return(mod)
}

