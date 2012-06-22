#require(auteur)
#require(diversitree)
#atr=as.integer(unlist(strsplit(toString(packageVersion("auteur")),".",fixed=TRUE))[3])
#div=as.integer(unlist(strsplit(toString(packageVersion("diversitree")),".",fixed=TRUE))[3])
#if(atr<507) {
#	stop("Update 'auteur': contact jonathan.eastman@gmail.com")
#}
#if(div<4){
# MUST BE DIVERSITREE 0.9-4 or greater due to bug fix in make.bm()
#	stop("Update 'diversitree': download and install from 'https://github.com/richfitz/diversitree'")
#}

#source("fitDiscrete.R")
#source("fitContinuous.R")

#source("sourcer.R")
### EXAMPLES ###
#require(auteur)


gen=FALSE
if(gen){
	require(auteur)
	phy=rescaleTree(ultrametricize.phylo(rtree(20), "max"),100)
	dis=as.integer(rTraitDisc(phy, matrix(c(-.1,.1,.1,.1,-.1,.1,.1,.1,-.1), nrow=3), states=c(1,2,3)))
	names(dis)=phy$tip.label
	con=rTraitCont(phy)
	m=matrix(sample(1:3,9, replace=TRUE), nrow=3) 
	mkn=make.mkn(phy,dis,k=length(unique(dis)))
	
	cat("\n\n\n*** 'm', 'phy','dis',and, 'con' are in memory ***\n")
	
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
	control=list(method=c("SANN","L-BFGS-B"), niter=100, FAIL=1e200, hessian=FALSE),
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
				   fd(phy, dat[,idx], model=model, transform=transform, bounds=bounds, control=control, ...)	   
				   })
		names(res)=nm
		return(res)
	} else {
		tmp=as.integer(dat[,1])
		names(tmp)=rownames(dat)
		dat=tmp
		if(!all(is.integer(dat))) stop("supply 'dat' as a vector (or matrix) of positive integers")

	}
	
	constrain=model
	model=match.arg(transform, c("none", "EB", "lambda", "kappa", "delta", "white"))

	# CONTROL OBJECT
	ct=list(method=c("SANN","L-BFGS-B"), niter=ifelse(model=="none", 50, 100), FAIL=1e200)
	if("method"%in%names(control)) control$method=match.arg(control$method, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), several.ok=TRUE)
	ct[names(control)]=control
	if(ct$niter<2) stop("'niter' must be equal to or greater than 2")
	
	lik=mkn.lik(phy, dat, constrain=constrain, transform=model, control=list(method="exp", root=ROOT.OBS), ...)
	if(model=="white") return(list(opt=lik))
	argn=unlist(argnames(lik))
	
	## CONSTRUCT BOUNDS ##
	mn=c(-5, -500, -500, -5, -500)
	mx=c(0.1, 0, 0, log(2.999999), log(1000))
	bnds=as.data.frame(cbind(mn, mx))
	bnds$typ=c("nat", "exp", "exp", "exp", "exp")
	rownames(bnds)=c("a", "lambda", "kappa", "delta", "trns")
	bnds$model=c("EB", "lambda", "kappa", "delta", "trns")
	
	parnm=ifelse(argn%in%rownames(bnds), argn, "trns")
	partp=bnds[parnm, "typ"]
	
	# User bounds
	if(length(bounds)>0){
		mm=match(names(bounds), rownames(bnds))
		if(any(is.na(mm)->zz)) mm[zz]=match(names(bounds)[zz], bnds$model)
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

	xx=function(p){
		if(par%in%c("a")){
			tmp=-lik(c(p[1], exp(p[-1])))
		} else {
			tmp=-lik(exp(p))
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
		typs=bnds[parnm, "typ"]
		
		# resolve method	
		if(length(argn)==1) {
			method="Brent" 
		} else {
			if(runif(1)<0.5){  ### FIXME ###
				method="subplex"
			} else {
				method=sample(ct$method,1)
			}
		}
		
		if(method=="subplex"){
			op<-out[[i]]<-try(suppressWarnings(subplex(par=start, fn=f, control=list(reltol = .Machine$double.eps^.25, parscale = rep(0.1, length(argn))), hessian=ct$hessian)),silent=TRUE)
		} else {
			op<-out[[i]]<-try(suppressWarnings(optim(par=start, fn=f, upper=max, lower=min, method=method, hessian=ct$hessian)),silent=TRUE)
		}
		
		if(!inherits(op,"try-error")){
			op$value=-op$value
			names(op)[names(op)=="value"]="lnL"
			names(op$par)=argn
			
			op$par=sapply(1:length(typs), function(x) if(typs[x]=="exp") return(exp(op$par[x])) else return(op$par[x]))
			mm[i, ]=c(op$par, op$lnL, op$convergence)
			mt[i]=method
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
		mod=mt[z]
	} else {
		mod=NA
		z=NA
		mm=c(rep(NA, length(argn)), -Inf)
		names(mm)=c(argn,"lnL")
	}
	k=length(argn)+1
	llx=which(names(mm)=="lnL")
	ll=mm[[llx]]
	mm=mm[-llx]
	
	## HESSIAN-based CI of par estimates
	if(ct$hessian){
		print(mm)
		print(partp)
		print(out[[z]]$hessian)
		hessian=out[[z]]$hessian
		CI=.bnd.hessian(hessian, mm, partp)
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
		allq=mm[argnames(lik)]
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
	
#	if(is.null(CI)){
		return(list(lik=lik, bnd=range[,c("mn", "mx")], res=res, opt=mm))	

#	} else {
#		return(list(lik=lik, bnd=range[,c("mn", "mx")], res=res, opt=mm, CI=CI))	

#	}
	
}

## WORKHORSE -- built from diversitree:::make.mkn  
# likelihood function creation
mkn.lik=function(
	phy, 
	dat, 
	constrain=c("ER","SYM","ARD","meristic"),
	transform=c("none", "EB","lambda", "kappa", "delta", "white"), 
	control=list(root=ROOT.OBS),
	...)
{
	
	# control object for make.mkn()
	ct = list(method="exp", root=ROOT.OBS)
	ct[names(control)]=control
	if(ct$method!="exp") stop(paste("method",sQuote(ct$method),"is not currently supported",sep=" "))
		
	# primary cache
	k<-nlevels(as.factor(dat))
    control <- check.control.mkn(ct, k)
    cache <- make.cache.mkn(phy, dat, k, strict=TRUE, control=ct)
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
	
	## KIND OF WORKING
	ll.mkn=function(cache, control) {
		k <- cache$info$k
		f.pars <- make.pars.mkn(k)
		f.pij <- make.pij.mkn(cache$info, control)
		idx.tip <- cache$idx.tip
		n.tip <- cache$n.tip
		n <- length(cache$len)
		map <- t(sapply(1:k, function(i) (1:k) + (i - 1) * k))
		idx.tip <- cbind(c(map[cache$states, ]), rep(seq_len(n.tip), k))
		children.C <- toC.int(t(cache$children))
		order.C <- toC.int(cache$order)
		
		.ll.mkn.exp=function(q, pars, intermediates=FALSE, preset = NULL) { # based on diversitree:::make.all.branches.mkn.exp
			if(is.null(argnames(FUN))) new=FUN() else new=FUN(q)
			
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
		if(is.null(argnames(FUN))){ # NO TRANSFORM
			ll=function(pars){
				qmat=f.pars(pars)
				ans=.ll.mkn.exp(q=NULL, pars=qmat, intermediates=FALSE)
				rootfunc.mkn(ans, qmat, ct$root, NULL, intermediates=FALSE)
			}
		} else {
			ll=function(pars){ # TREE TRANSFORM
				qmat=f.pars(pars[-1])
				ans=.ll.mkn.exp(q=pars[1], pars=qmat, intermediates=FALSE)
				rootfunc.mkn(ans, qmat, ct$root, NULL, intermediates=FALSE)
			}
			
		}
		class(ll) <- c("mkn", "dtlik", "function")
		attr(ll,"argnames") <- c(argnames(FUN), cache$info$argnames)
		return(ll)
	}
	
	lik=ll.mkn(cache, control)
	
	## CONSTRAINTS
	if(!all(constrain=="ARD")){
		if(is.character(constrain)){
			cc=match.arg(constrain, c("ER","SYM","ARD","meristic"))
			lik=constrain.k(lik, model=cc, ...)
		} else {
			if(is.matrix(constrain)){
				if(ncol(constrain)==max(dat)){
					lik=constrain.m(lik, m=constrain)
				}
			} else {
				stop("'constrain' must be supplied as a dummy matrix representing constraints on transition classes")
			}
		}
	}
	lik
}

