.treedata <-
function(phy, data)
{
	
#	if(is.null(data.names)) 
#	{
	if(is.vector(data)){
		data.names<-names(data)
	} else {
		data.names<-rownames(data)
	}
#	}
	t<-phy$tip.label
	r1<-t[is.na(match(t,data.names))]
	r2<-data.names[is.na(match(data.names,t))]
	
	r<-list(sort(r1), sort(r2))
	
	names(r)<-cbind("tree_not_data", "data_not_tree")
	if(length(r1)==0 && length(r2)==0) return("OK")
	else return(r)
}

cherries <- function(phy){
    cache=.cache.descendants(phy)
	N=Ntip(phy)
	nds=which(tabulate(phy$edge[phy$edge[,2]<=N,1])==2)
    if(length(nds)){
        mat=matrix(NA, nrow=length(nds), ncol=2)
        for(i in 1:length(nds)) mat[i, ]=phy$tip.label[cache$tips[[nds[i]]]]
        rownames(mat)=nds
    } else {
        mat=NULL
    }
    return(mat)
}


# Treedata is a function internal to GEIGER
# that makes sure that the names of the taxa 
# in the tree and data file match, and prunes 
# things accordingly; it returns a list with 
# two elements: phy and data

tips <- function(phy, node)
{
	if(node<=Ntip(phy)) return(phy$tip.label[node])
	dd=.get.descendants.of.node(node, phy, tips=TRUE)
	phy$tip.label[dd]
}

span.phylo=function(phy){
    desc=.cache.descendants(phy)
    N=Ntip(phy)
    labs=c(phy$tip.label, phy$node.label)
    
    ff=function(node){
        if(length(node)>1) stop("Supply 'node' as a single value")
        if(!is.numeric(node)) node=which(labs==node)
        
        if(node<=N) return(NULL)
        dd=desc$fdesc[[node]]
        labs[sapply(dd, function(x) desc$tips[[x]][1])]
    }
    ff
    
}

# data.names is optional, and will replace the names or rownames
# of data when matching data to the tree

# if sort is T, data will have rows in the same order
# as the taxon names in phy$tip.label

treedata<-function(phy, data, sort=FALSE, warnings=TRUE)
{
	
	dm=length(dim(data))
	
	if(is.vector(data)) {
		data<-as.matrix(data)
	}
	if(is.factor(data)) {
		data<-as.matrix(data)
	}
	if(is.array(data) & length(dim(data))==1) {
		data<-as.matrix(data)
	}
	
#	if(is.null(data.names)) {
	if(is.null(rownames(data))) {
		stop("names for 'data' must be supplied")
#JME				data.names<-phy$tip.label
#JME				if(warnings)
#JME					cat("Warning: no tip labels, order assumed to be the same as in the tree\n")
	} else {
		data.names<-rownames(data)
	}
#	}
	nc<-.treedata(phy, data)
	if(is.na(nc[[1]][1]) | nc[[1]][1]!="OK") {
		if(length(nc[[1]]!=0)) {
			phy=drop.tip(phy, as.character(nc[[1]]))
			if(warnings) {
				warning(paste("The following tips were not found in 'data' and were dropped from 'phy':\n\t",
							  paste(nc[[1]], collapse="\n\t"), sep=""))
#JME			print(nc[[1]])
#JME			cat("\n")
			}
		}
		
		if(length(nc[[2]]!=0)) {
			m<-match(data.names, nc[[2]])
			data=as.matrix(data[is.na(m),])
			data.names<-data.names[is.na(m)]
			if(warnings) {
				warning(paste("The following tips were not found in 'phy' and were dropped from 'data':\n\t",
							  paste(nc[[2]], collapse="\n\t"), sep=""))
#JME			print(nc[[2]])
#JME			cat("\n")
			}
		}
 	}
	order<-match(data.names, phy$tip.label)	
	
	rownames(data)<-phy$tip.label[order]
	
	if(sort) {
		
    	index <- match(phy$tip.label, rownames(data))
   		data <- as.matrix(data[index,])
	} 
	if(dm==2){
		data <- as.matrix(data)
	}
	
	phy$node.label=NULL
	
	return(list(phy=phy, data=data))
}

## GENERIC 
transform=function(x, ...) UseMethod("transform")


## GENERIC
print.transformer=function (x, printlen = 3, ...) 
{
    cat("function for tree transformation\n")
    cat("\targument names:", paste(argn(x), collapse = ", "))		 
	cat("\n\n")
	f=x
	attributes(f)=NULL
	cat("definition:\n")
	print(f)
}



## GENERIC
print.bm=function (x, printlen = 3, ...) 
{
    cat("likelihood function for univariate continuous trait evolution\n")
    cat("\targument names:", paste(argn(x), collapse = ", "))		 
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
    cat("likelihood function for univariate discrete trait evolution\n")
    cat("\targument names:", paste(argn(x), collapse = ", "))	
    if(!is.null(al<-attr(x, "levels"))) {
        fmt=.format.levels.print(length(al))
        cat("\n\n\tmapping\n\t\t", paste(sprintf(fmt, 1:length(al)), al, collapse="\n\t\t", sep=": "), sep="")
    }
	cat("\n\n")
	f=x
	attributes(f)=NULL
	cat("definition:\n")
	print(f)
}

## GENERIC
argn.default=function(x, ...){
	attr(x, "argn")
}

## GENERIC
argn.mkn=function(x, ...){
	attr(x, "argn")
}

## HESSIAN COMPUTATION of CONFIDENCE INTERVALS
.bnd.hessian=function(m, p, s, prob=0.05){
	prob=min(c(1-prob, prob))
	dm=unique(dim(m))
	if(length(dm)!=1) {warning("FAILURE in Hessian CI computation: 'm' must be a square matrix"); return(NA)}
	if(dm!=length(p) | dm!=length(s)) {warning("FAILURE in Hessian CI computation: 'p' and 's' must be of equal length to the dimensions of 'm'"); return(NA)}
	if(!all(s%in%c("exp","nat"))) {warning("FAILURE in Hessian CI computation: 's' must indicate the space of 'p': expecting either 'exp' or 'nat' as elements"); return(NA)}
	
	qq=qnorm(1-(prob/2))
	dd=try(sqrt(diag(solve(m))),silent=TRUE)
	if(inherits(dd, "try-error")){
		warning(paste("ERROR in inversion of Hessian matrix:", gsub(">", "", gsub("<", "", gsub("\n", "", toString(attributes(dd)$condition))))))
		return(NA)
	}
	bb=qq*dd
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

.repars=function(pars, expected){
    if(!length(pars)==length(expected)) stop(paste("The following 'pars' are expected:\n\t", paste(expected, collapse="\n\t", sep=""), sep=""))
	if(all(!is.null(nm<-names(pars)))){
		if(!all(nm%in%expected)) stop(paste("The following 'pars' are unexpected:\n\t", paste(nm[!nm%in%expected], collapse="\n\t", sep=""), sep=""))
		if(length(unique(nm))!=length(expected)) stop(paste("The following 'pars' are expected:\n\t", paste(expected, collapse="\n\t", sep=""), sep="")) 
		mm=match(expected, nm)
		return(pars[mm])
	} else {
		return(pars)
	}
}

## EXTENSION of diversitree:::constrain()
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

## EXTENSION of diversitree:::constrain()
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
.ou.smartstart=function(dat, bounds){
	vv=var(dat)
    xb=max(bounds)
    nb=min(bounds)
    atry=seq(-8,4,by=2)
    s=sample(1:length(atry),1)
    if(s==1) {
        aa=nb
    } else if(s==length(atry)) {
        aa=xb
    } else {
        aa=vv*2*exp(atry[s])
        if(aa>xb) aa=xb
        if(aa<nb) aa=nb
    }
    if(is.na(aa)) aa=0.1
    aa
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


# compute path length from root to tip
.paths.cache=function(cache){
	
    if(is.null(cache$ordering) || cache$ordering!="postorder"){
        stop("'cache' should be postordered")
    }

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
	phy <- reorder(phy, "postorder")

	n <- length(phy$tip.label)
	pp <- .cache.descendants(phy)$adesc[-c(1:n)]
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



## DEFUNCT
.DEFUNCT_reorder.cache.pruningwise=function(cache){
	if(cache$ordering!="pruningwise") {
		nb.node <- cache$n.node
		if (nb.node != 1) {
			nb.tip <- cache$n.tip
			nb.edge <- nrow(cache$edge)
			ord=.C("neworder_pruningwise", as.integer(nb.tip), as.integer(nb.node), 
				   as.integer(cache$edge[, 1]), as.integer(cache$edge[, 2]), as.integer(nb.edge), 
				   integer(nb.edge), PACKAGE = "geiger")[[6]]
			cache$edge=cache$edge[ord,]
			cache$edge.length=cache$edge.length[ord]
		}	
		cache$ordering="pruningwise"
	}
	cache
}

## DEFUNCT
.DEFUNCT_reorder.cache.cladewise=function(cache){
	if(cache$ordering!="cladewise") {
		nb.node <- cache$n.node
		if (nb.node != 1) {
			nb.tip <- cache$n.tip
			nb.edge <- nrow(cache$edge)
			ord=.C("neworder_phylo", as.integer(nb.tip), as.integer(cache$edge[, 1]), 
				   as.integer(cache$edge[, 2]), as.integer(nb.edge), integer(nb.edge), as.integer(1), DUP=FALSE, NAOK=TRUE, PACKAGE="geiger")[[5]]
			cache$edge=cache$edge[ord,]
			cache$edge.length=cache$edge.length[ord]
		}	
		cache$ordering="cladewise"
	}
	cache
}

.heights.cache=function (cache) 
{
    if(is.null(cache$ordering) || cache$ordering!="postorder"){
        stop("'cache' should be postordered")
    }
        
	n <- cache$n.tip
	n.node <- cache$n.node
	xx <- numeric(n + n.node)
	for (i in nrow(cache$edge):1) xx[cache$edge[i, 2]] <- xx[cache$edge[i, 1]] + cache$edge.length[i]
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

## tree transformation
transform.phylo=function(x, model=c("OU", "EB", "nrate", "lrate", "trend", "lambda", "kappa", "delta", "white", "depth"), ...){
	
	phy=x
	
	model=match.arg(model, c("OU", "EB", "nrate", "lrate", "trend", "lambda", "kappa", "delta", "white", "depth"))
	
	if(!"phylo"%in%class(phy)) stop("supply 'phy' as a 'phylo' object")
	
	FUN=switch(model, 
			   OU=.ou.phylo(phy),
			   EB=.eb.phylo(phy),
               nrate=.nrate.phylo(phy),
               lrate=.lrate.phylo(phy),
			   trend=.trend.phylo(phy),
			   lambda=.lambda.phylo(phy),
			   kappa=.kappa.phylo(phy),
			   delta=.delta.phylo(phy),
			   white=.white.phylo(phy),
			   depth=.depth.phylo(phy)
			   )
	class(FUN)=c("transformer", "function")
	dots=list(...)
	if(length(dots)==length(argn(FUN))) return(FUN(unlist(dots)))
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
.nrate.phylo=function(phy){
	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t=Tmax-ht$end
	ht$e=ht$start-ht$end
	ht$a=ht$t-ht$e
	ht$rS=ht$a/Tmax
	ht$rE=ht$t/Tmax
	
	dd=phy$edge[,2]
	
	relscale.brlen=function(start, end, len, dat){
		ss=start<=dat[,"time"]
		strt=min(which(ss))
		
		ee=dat[,"time"]<end
		etrt=max(which(ee))+1
		
		bl=numeric()
        fragment=numeric()
		marker=start
		for(i in strt:etrt){
			fragment=c(fragment, (nm<-(min(c(end, dat[i, "time"]))))-marker)
			bl=c(bl,dat[i, "rate"])
			marker=nm
		}
        fragment=fragment/(sum(fragment))
        sclbrlen=numeric()
        for(i in 1:length(bl)) sclbrlen=c(sclbrlen, len*fragment[i]*bl[i])
		sc=structure(as.numeric(sclbrlen), names=strt:etrt)
		return(sc)
	}
	
	
	z=function(time, rate, rescale=TRUE){
		if(any(time>1) | any(time<0)) stop("supply 'time' as a vector of relative time:\n\tvalues should be in the range 0 (root) to 1 (present)")
		if(any(rate<0)) stop("'rate' must consist of positive values")
		if(length(time)!=length(rate)) stop("'time' and 'rate' must be of equal length")
		ordx=order(time)
		time=time[ordx]
		rate=rate[ordx]
		dat=cbind(time=c(0,time, 1), rate=(c(1, 1, rate)))
		rs=sapply(dd, function(x) as.numeric(sum(relscale.brlen(ht$rS[x], ht$rE[x], ht$e[x], dat))))
		phy$edge.length=rs
        if(rescale) phy=transform.phylo(phy, "depth", Tmax)
		phy
	}
	attr(z, "argn")=c("time", "rate")
	return(z)
}

.lrate.phylo=function(phy){
    N=Ntip(phy)
    n=Nnode(phy)
    cache=.cache.tree(phy)
    cache$phy=phy
    vv=rep(1, N+n-1)

	z=function(node, rate){
        
        
        shifts=c(sort(node[node>N]), node[node<=N])
        mm=match(shifts, node)
        rates=rate[mm]
        
        for(i in 1:length(shifts)) vv=.assigndescendants(vv, shifts[i], rates[i], exclude=shifts, cache=cache)
        phy$edge.length=vv*phy$edge.length
        phy
    }
    
	attr(z, "argn")=c("node", "rate")
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
	attr(z,"argn")="lambda"
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
	attr(z,"argn")="lambda"
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
	attr(z,"argn")="delta"
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
	attr(z,"argn")="delta"
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
        scl=sapply(1:nrow(ht), function(idx){
            if(idx==N+1) return(NA)
            if(ht$br[idx]>0 & ht$er[idx]>0) {
                return((ht$br[idx]+ht$er[idx])/2)
            } else if(ht$br[idx]<0 & ht$er[idx]<0) {
                return(0)
            } else {
                si=-1/slope
                return(ht$br[idx]*(si-ht$head[idx])/(2*(ht$tail[idx]-ht$head[idx])))
            }
        })
		cache$len=cache$len*scl
		cache
	}
	attr(z,"argn")="slope"
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
        # begin (age): head
        # end: tail
		ht$br=1+ht$head*slope
		ht$er=1+ht$tail*slope
        scl=sapply(1:nrow(ht), function(idx){
            if(idx==N+1) return(NA)
            if(ht$br[idx]>0 & ht$er[idx]>0) {
                return((ht$br[idx]+ht$er[idx])/2)
            } else if(ht$br[idx]<0 & ht$er[idx]<0) {
                return(0)
            } else {
                si=-1/slope
                return(ht$br[idx]*(si-ht$head[idx])/(2*(ht$tail[idx]-ht$head[idx])))
            }
        })
        
		phy$edge.length=phy$edge.length*scl[phy$edge[,2]]
		phy
	}
	attr(z,"argn")="slope"
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
	attr(z,"argn")="alpha"
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
	attr(z,"argn")="alpha"
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
	attr(z,"argn")="a"
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
	attr(z,"argn")="a"
	return(z)
}



# tree transformation
.kappa.cache=function(cache){
	
	z=function(kappa){
		if(kappa<0) stop("'kappa' must be positive valued")
		
		cache$len=cache$len^kappa
		cache
	}
	attr(z,"argn")="kappa"
	return(z)
}

# tree transformation
.kappa.phylo=function(phy){
	
	z=function(kappa){
		if(kappa<0) stop("'kappa' must be positive valued")
		
		phy$edge.length=phy$edge.length^kappa
		phy
	}
	attr(z,"argn")="kappa"
	return(z)
}

# tree transformation
.depth.phylo=function(phy){
	orig=max(heights.phylo(phy))
	z=function(depth){
		phy$edge.length <- (phy$edge.length/orig) * depth
		if(!is.null(phy$root.edge)) phy$root.edge=(phy$root.edge/orig) * depth
		phy
	}
	attr(z,"argn")="depth"
	z
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

drop.extinct<-function (phy, tol = .Machine$double.eps^0.5) 
{
    if (class(phy) != "phylo") 
	stop("object \"phy\" is not of class \"phylo\"")
    phy2 <- phy
    phy <- new2old.phylo(phy)
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    xx <- as.numeric(rep(NA, nb.tip + nb.node))
    names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
    xx["-1"] <- 0
    for (i in 2:length(xx)) {
        nod <- names(xx[i])
        ind <- which(as.numeric(phy$edge[, 2]) == nod)
        base <- phy$edge[ind, 1]
        xx[i] <- xx[base] + phy$edge.length[ind]
    }
    depth <- max(xx)
    offset <- depth - xx[names(xx) > 0]
    drops <- phy$tip.label[offset > tol]
    if (length(drops) >= (nb.tip - 1)) 
	return(NULL)
    if(length(drops)==0)
	return(phy2)    
    res <- drop.tip(phy2, drops)
    res
}

drop.random<-function (phy, n) 
{
    if (class(phy) != "phylo") 
	stop("object \"phy\" is not of class \"phylo\"")
    nb.tip <- length(phy$tip.label)
    if (n > nb.tip) 
	return(NULL)
    cut <- sample(1:nb.tip, n)
    r <- drop.tip(phy, cut)
    return(r)
}

