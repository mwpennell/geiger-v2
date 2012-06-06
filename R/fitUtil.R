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


### EXAMPLES ###

#phy=rescaleTree(ultrametricize.phylo(rtree(2000), "max"),100)
#dis=as.integer(rTraitDisc(phy, matrix(c(-.1,.1,.1,.1,-.1,.1,.1,.1,-.1), nrow=3), states=c(1,2,3)))
#names(dis)=phy$tip.label
#con=rTraitCont(phy)
#m=matrix(sample(1:3,9, replace=TRUE), nrow=3) 
#mkn=make.mkn(phy,dis,k=length(unique(dis)))

#cat("\n\n\n*** 'm', 'phy','dis',and'con' are in memory ***\n")

###############


.bnd.hessian=function(m, p, s){
	dm=unique(dim(m))
	if(length(dm)!=1) {warning("FAILURE in Hessian CI computation: 'm' must be a square matrix"); return(NA)}
	if(dm!=length(p) | dm!=length(s)) {warning("FAILURE in Hessian CI computation: 'p' and 's' must be of equal length to the dimensions of 'm'"); return(NA)}
	if(!all(s%in%c("exp","nat"))) {warning("FAILURE in Hessian CI computation: 's' must indicate the space of 'p': expecting either 'exp' or 'nat' as elements"); return(NA)}

	qq=qnorm(0.975)
	dd=try(sqrt(diag(solve(m))),silent=TRUE)
	if(inherits(dd, "try-error")){
		warning(paste("ERROR in inversion of Hessian matrix:", gsub(">", "", gsub("<", "", gsub("\n", "", toString(attributes(dd)$condition))))))
		return(NA)
	}
	bb=qq*dd
#	bb=ifelse(s=="exp", exp(qq)*dd, qq*dd)
	res=sapply(1:length(s), function(idx) {
		if(s[idx]=="exp"){
		   
		   return(c(lb=exp(log(p[idx])-bb[idx]), ub=exp(log(p[idx])+bb[idx])))
		} else {
		   return(c(lb=p[idx]-bb[idx], ub=p[idx]+bb[idx]))
		}
	})
	rownames(res)=c("lb", "ub")
	colnames(res)=names(p)
	res
	
}


## MKN CONSTRAINT FUNCTIONS
.er.matrix=function(m){
	m[]=1
	diag(m)=0
	m
}

.sym.matrix=function(m){
	m[]=1:length(m)
	ww=which(upper.tri(m),arr.ind=TRUE)
	m[ww]=m[ww[,2:1]]
	diag(m)=0
	m
}

.ard.matrix=function(m){
	m[]=1:length(m)
	diag(m)=0
	m
}

.meristic.matrix=function(m, symmetric=TRUE){
	m[]=0
	k=nrow(m)
	idx=rbind(cbind(1:(k-1), 2:k), cbind(2:k, 1:(k-1)))
	if(symmetric) base=.sym.matrix(m) else base=.ard.matrix(m)
	m[idx]=base[idx]
	m
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
	return(constrain.m(f, mm))
}


## create constrained likelihood function by supplying matrix
constrain.m=function(f, m){
	# if entry is 0, assumed to constrain rate to 0
	k=ncol(m)
	if(nrow(m)!=k) stop("supply 'm' as a square matrix")	
	idx <- cbind(rep(1:k, each = k - 1), unlist(lapply(1:k, function(i) (1:k)[-i])))
	npar <- k * (k - 1)
	par=cbind(idx, m[idx], duplicated(m[idx]))
	colnames(par)=c("row","col","parm","dup")
	upar=unique(par[,"parm"])
	sdig=nchar(k)
	sfor=paste("%0",sdig,"d",sep="")
	spar=apply(par[,c("row","col")], 1, function(x) paste("q", paste(sapply(x, function(y) sprintf(sfor, y)), collapse=""),sep=""))
	res=unlist(sapply(1:nrow(par), function(idx){
				  curspar=spar[idx]
				  if(par[idx,"parm"]==0){
					return(paste(curspar, 0, sep="~") )
				  } else if(par[idx,"dup"]==1){
					mm=min(which(par[,"parm"]==par[idx,"parm"]))
					return(paste(curspar, spar[mm], sep="~") )
				  } else {
					return(NULL)
				  }  	
				  })) 
	return(constrain(f, formulae=res))
}

# (smart) starting pt for optimization
.ou.smartstart=function(sigsq, var){
	alpha=(2*sigsq)/var
	alpha
}

# (smart) starting pt for optimization
.bm.smartstart=function(phy, dat){
	ss=mean(pic(dat, phy)^2)
	ss
}


.aic=function(v, n){
# v: has object 'lnL' and 'k'
	v$aic <- 2 * v$k - 2 * v$lnL
    v$aicc <- 2 * v$k * (n - 1)/(n - v$k - 2) - 2 * v$lnL
	v
}

## GENERIC
print.transformer=function (x, printlen = 3, ...) 
{
    cat("function for tree transformation\n")
    cat("\targnames:", paste(argnames(x), collapse = ", "))		 
	cat("\n\n")
	f=x
	attributes(f)=NULL
	cat("definition:\n")
	print(f)
}



## GENERIC
print.bm=function (x, printlen = 3, ...) 
{
    cat("likelihood function for generalized branching diffusion\n")
    cat("\targnames:", paste(argnames(x), collapse = ", "))		 
	cat("\n\n")
	f=x
	attributes(f)=NULL
	cat("definition:\n")
	print(f)
}

## GENERIC
print.mkn=
function (x, printlen = 3, ...) 
{
    cat("likelihood function for Mk(n)\n")
    cat("\targnames:", paste(argnames(x), collapse = ", "))		 
	cat("\n\n")
	f=x
	attributes(f)=NULL
	cat("definition:\n")
	print(f)
}



## GENERIC
argnames.default=function(x, ...){
	attr(x, "argnames")
}

## GENERIC
argnames.mkn=function(x, ...){
	attr(x, "argnames")
}


# compute path length from root to tip
.paths.cache=function(cache){
	cache=.reorder.cache.pruningwise(cache)
	n <- cache$n.tip
	pp <- cache$adesc[-c(1:n)]
	e1 <- cache$edge[, 1]
	e2 <- cache$edge[, 2]
	EL <- cache$edge.length
	xx <- numeric(n + cache$n.node)
	for (i in length(e1):1) {
		var.cur.node <- xx[e1[i]]
		xx[e2[i]] <- var.cur.node + EL[i]
		j <- i - 1L
		while (e1[j] == e1[i] && j > 0) {
			left <- if (e2[j] > n) 
			pp[[e2[j] - n]]
			else e2[j]
			right <- if (e2[i] > n) 
			pp[[e2[i] - n]]
			else e2[i]
			j <- j - 1L
		}
	}
	xx[1:n]
}

# compute path length from root to tip
.paths.phylo=function(phy){
## from vcv.phylo()
	n <- length(phy$tip.label)
	pp <- auteur:::.compile_descendants(phy)$adesc[-c(1:n)]
	phy <- reorder(phy, "pruningwise")
	e1 <- phy$edge[, 1]
	e2 <- phy$edge[, 2]
	EL <- phy$edge.length
	xx <- numeric(n + phy$Nnode)
	for (i in length(e1):1) {
		var.cur.node <- xx[e1[i]]
		xx[e2[i]] <- var.cur.node + EL[i]
		j <- i - 1L
		while (e1[j] == e1[i] && j > 0) {
			left <- if (e2[j] > n) 
			pp[[e2[j] - n]]
			else e2[j]
			right <- if (e2[i] > n) 
			pp[[e2[i] - n]]
			else e2[i]
			j <- j - 1L
		}
	}
	xx[1:n]
}




## tree transformation
transform.phylo=function(phy, model=c("OU", "EB", "trend", "lambda", "kappa", "delta", "white")){
	
	require(auteur)
	
	model=match.arg(model, c("OU", "EB", "trend", "lambda", "kappa", "delta", "white"))
	
	if(!"phylo"%in%class(phy)) stop("supply 'phy' as a 'phylo' object")
	
	FUN=switch(model, 
			   OU=.ou.phylo(phy),
			   EB=.eb.phylo(phy),
			   trend=.trend.phylo(phy),
			   lambda=.lambda.phylo(phy),
			   kappa=.kappa.phylo(phy),
			   delta=.delta.phylo(phy),
			   white=.white.phylo(phy)
			   )
	class(FUN)=c("transformer", "function")
	return(FUN)
	
}

.reorder.cache.pruningwise=function(cache){
	if(cache$ordering!="pruningwise") {
		nb.node <- cache$n.node
		if (nb.node != 1) {
			nb.tip <- cache$n.tip
			nb.edge <- dim(cache$edge)[1]
			ord=.C("neworder_pruningwise", as.integer(nb.tip), as.integer(nb.node), 
				   as.integer(cache$edge[, 1]), as.integer(cache$edge[, 2]), as.integer(nb.edge), 
				   integer(nb.edge), PACKAGE = "ape")[[6]]
			cache$edge=cache$edge[ord,]
			cache$edge.length=cache$edge.length[ord]
		}	
		cache$ordering="pruningwise"
	}
	cache
}

.reorder.cache.cladewise=function(cache){
	if(cache$ordering!="cladewise") {
		nb.node <- cache$n.node
		if (nb.node != 1) {
			nb.tip <- cache$n.tip
			nb.edge <- dim(cache$edge)[1]
			ord=.C("neworder_cladewise", as.integer(nb.tip), as.integer(x$edge[, 1]), 
				   as.integer(x$edge[, 2]), as.integer(nb.edge), integer(nb.edge), PACKAGE = "ape")[[5]]
			cache$edge=cache$edge[ord,]
			cache$edge.length=cache$edge.length[ord]
		}	
		cache$ordering="cladewise"
	}
	cache
}

.heights.cache<-
function (cache) 
{
	cache=.reorder.cache.cladewise(cache)
	n <- cache$n.tip
	n.node <- cache$n.node
	xx <- numeric(n + n.node)
	for (i in 1:nrow(cache$edge)) xx[cache$edge[i, 2]] <- xx[cache$edge[i, 1]] + cache$edge.length[i]
	root = ifelse(is.null(cache$root.edge), 0, cache$root.edge)
	depth = max(xx)
	tt = depth - xx
	idx = 1:length(tt)
	dd = cache$edge.length[idx]
	mm = match(1:length(tt), c(cache$edge[, 2], n + 1))
	dd = c(cache$edge.length, root)[mm]
	ss = tt + dd
	res = cbind(ss, tt)
	rownames(res) = idx
	colnames(res) = c("start", "end")
	res = data.frame(res)
	res	
	
}

## REPLACE in AUTEUR (move to heights.phylo)
.heights.phylo=function(phy){
	phy <- reorder(phy, order="cladewise")
	n <- length(phy$tip.label)
	n.node <- phy$Nnode
	xx <- numeric(n + n.node)
	for (i in 1:nrow(phy$edge)) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	root = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
	depth = max(xx)
	tt = depth - xx
	idx = 1:length(tt)
	dd = phy$edge.length[idx]
	mm = match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
	dd = c(phy$edge.length, root)[mm]
	ss = tt + dd
	res = cbind(ss, tt)
	rownames(res) = idx
	colnames(res) = c("start", "end")
	res = data.frame(res)
	res	
}

transform.phylo=function(phy, model=c("OU", "EB", "trend", "lambda", "kappa", "delta", "white")){
	model=match.arg(model, c("OU", "EB", "trend", "lambda", "kappa", "delta", "white"))
	
	if(!"phylo"%in%class(phy)) stop("supply 'phy' as a 'phylo' object")
	
	FUN=switch(model, 
			   OU=.ou.phylo(phy),
			   EB=.eb.phylo(phy),
			   trend=.trend.phylo(phy),
			   lambda=.lambda.phylo(phy),
			   kappa=.kappa.phylo(phy),
			   delta=.delta.phylo(phy),
			   white=.white.phylo(phy)
			   )
	
	return(FUN)
}



# tree transformation
.white.cache=function(cache){
	N=cache$n.tip
	cache$len[]=0
	cache$len[1:N]=1
	
	z=function(){
		cache
	}
	return(z)
}

# tree transformation
.white.phylo=function(phy){
	N=Ntip(phy)
	phy$edge.length[]=0
	phy$edge.length[phy$edge[,2]<=N]=1
	
	z=function(){
		phy
	}
	return(z)
}



# tree transformation
.null.cache=function(cache) {
	z=function(){
		cache
	}
	return(z)
}

# tree transformation
.null.phylo=function(phy) {
	z=function(){
		phy
	}
	return(z)
}



# tree transformation
.lambda.cache=function(cache){
	N=cache$n.tip
	paths=.paths.cache(cache)
	
	z=function(lambda){
		if(lambda<0) stop("'lambda' must be positive valued")
		
		bl=cache$len*lambda
		bl[1:N]=bl[1:N]+(paths-(paths*lambda))
		cache$len=bl
		if(any(cache$len<0,na.rm=TRUE)){
			warning("negative branch lengths encountered:\n\tlambda may be too large")
		}
		cache
	}
	attr(z,"argnames")="lambda"
	return(z)
}

# tree transformation
.lambda.phylo=function(phy){
	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:N, phy$edge[,2])
	ht$e=ht$start-ht$end
	paths=.paths.phylo(phy)
	
	z=function(lambda){
		if(lambda<0) stop("'lambda' must be positive valued")
		
		bl=phy$edge.length*lambda
		bl[mm]=bl[mm]+(paths-(paths*lambda))
		phy$edge.length=bl
		if(any(phy$edge.length<0)){
			warning("negative branch lengths encountered:\n\tlambda may be too large")
		}
		phy
	}
	attr(z,"argnames")="lambda"
	return(z)
}



# tree transformation
.delta.cache=function(cache){
	ht=.heights.cache(cache)
	N=cache$n.tip
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), cache$edge[,2])
	ht$t=Tmax-ht$end
	ht$e=ht$start-ht$end
	ht$a=ht$t-ht$e
	
	z=function(delta, rescale=TRUE){
		if(delta<0) stop("'delta' must be positive valued")
		bl=(ht$a+ht$e)^delta-ht$a^delta
		cache$len=bl
		
		if(rescale){
			scl=Tmax^delta
			cache$len=(cache$len/scl)*Tmax
		}
		cache
	}
	attr(z,"argnames")="delta"
	return(z)
}

# tree transformation
.delta.phylo=function(phy){
	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t=Tmax-ht$end
	ht$e=ht$start-ht$end
	ht$a=ht$t-ht$e
	
	z=function(delta, rescale=TRUE){
		if(delta<0) stop("'delta' must be positive valued")
		bl=(ht$a+ht$e)^delta-ht$a^delta
		phy$edge.length=bl[phy$edge[,2]]
		
		if(rescale){
			scl=Tmax^delta
			phy$edge.length=(phy$edge.length/scl)*Tmax
		}
		phy
	}
	attr(z,"argnames")="delta"
	return(z)
}



# tree transformation
.trend.cache=function(cache){
	
	ht=.heights.cache(cache)
	N=cache$n.tip
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), cache$edge[,2])
	ht$head=Tmax-ht$end[cache$edge[mm,1]] # age
	ht$tail=ht$head+(ht$start-ht$end)
	
	
	z=function(slope){
		ht$br=1+ht$head*slope
		ht$er=1+ht$tail*slope
		bl=ifelse(ht$br>0 & ht$er>0, 
				  (ht$br+ht$er)/2, 
				  ifelse(ht$br<0 & ht$er<0, 
						 0, 
						 ht$br*((-1/slope)-ht$head)/(2*(ht$tail-ht$head))))
		cache$len=cache$len*bl
		cache
	}
	attr(z,"argnames")="slope"
	return(z)
}

# tree transformation
.trend.phylo=function(phy){
	
	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$head=Tmax-ht$end[phy$edge[mm,1]] # age
	ht$tail=ht$head+(ht$start-ht$end)
	
	
	z=function(slope){
		ht$br=1+ht$head*slope
		ht$er=1+ht$tail*slope
		bl=ifelse(ht$br>0 & ht$er>0, 
				  (ht$br+ht$er)/2, 
				  ifelse(ht$br<0 & ht$er<0, 
						 0, 
						 ht$br*((-1/slope)-ht$head)/(2*(ht$tail-ht$head))))
		phy$edge.length=phy$edge.length*bl[phy$edge[,2]]
		phy
	}
	attr(z,"argnames")="slope"
	return(z)
}



# tree transformation
.ou.cache=function(cache){
	ht=.heights.cache(cache)
	N=cache$n.tip
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), cache$edge[,2])
	ht$t1=Tmax-ht$end[cache$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1
	z=function(alpha){
		if(alpha<0) stop("'alpha' must be positive valued")
		bl=(1/(2*alpha))*exp(-2*alpha * (Tmax-ht$t2)) * (1 - exp(-2 * alpha * ht$t2)) - (1/(2*alpha))*exp(-2*alpha * (Tmax-ht$t1)) * (1 - exp(-2 * alpha * ht$t1))
		cache$len=bl
		cache
	}
	attr(z,"argnames")="alpha"
	return(z)
}

# tree transformation
.ou.phylo=function(phy){
	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t1=Tmax-ht$end[phy$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1
	z=function(alpha){
		if(alpha<0) stop("'alpha' must be positive valued")
		bl=(1/(2*alpha))*exp(-2*alpha * (Tmax-ht$t2)) * (1 - exp(-2 * alpha * ht$t2)) - (1/(2*alpha))*exp(-2*alpha * (Tmax-ht$t1)) * (1 - exp(-2 * alpha * ht$t1))
		phy$edge.length=bl[phy$edge[,2]]
		phy
	}
	attr(z,"argnames")="alpha"
	return(z)
}



# tree transformation
.eb.cache=function(cache){
	
	ht=.heights.cache(cache)
	N=cache$n.tip
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), cache$edge[,2])
	ht$t1=Tmax-ht$end[cache$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1
	
	z=function(a){
		if(a==0) return(cache)
		bl = (exp(a*ht$t2)-exp(a*ht$t1))/(a)
		cache$len=bl
		cache
	}
	attr(z,"argnames")="a"
	return(z)
}

# tree transformation
.eb.phylo=function(phy){
	
	
	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t1=Tmax-ht$end[phy$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1
	
	z=function(a){
		if(a==0) return(phy)
		bl = (exp(a*ht$t2)-exp(a*ht$t1))/(a)
		phy$edge.length=bl[phy$edge[,2]]
		phy
	}
	attr(z,"argnames")="a"
	return(z)
}



# tree transformation
.kappa.cache=function(cache){
	
	z=function(kappa){
		if(kappa<0) stop("'kappa' must be positive valued")
		
		cache$len=cache$len^kappa
		cache
	}
	attr(z,"argnames")="kappa"
	return(z)
}

# tree transformation
.kappa.phylo=function(phy){
	
	z=function(kappa){
		if(kappa<0) stop("'kappa' must be positive valued")
		
		phy$edge.length=phy$edge.length^kappa
		phy
	}
	attr(z,"argnames")="kappa"
	return(z)
}


white.mkn=function(dat){
	tt=table(dat)
	n=sum(tt)
	p=tt/n
	ll=sum(log(p^tt))
	
	k=length(tt)-1
	opt=list(lnL=ll, method="MLE", k=k)
	opt=.aic(opt, n)
	opt
}

