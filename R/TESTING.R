.geigerwarn=function(...) warning("the called function is currently in development and is not fully vetted", ...)

.mcmc=function(lik, prior=list(), start=NULL, proposal=NULL, control=list(n=1e4, s=100, w=1)){
#mcmc(lik, start=c(sigsq=1, SE=1, z0=4), prior=list(sigsq=function(x) dexp(x, 1/1000, log=TRUE), SE=function(x) dexp(x, 1/1000, log=TRUE), z0=function(x) dnorm(x, mean=mean(dat), sd=100, log=TRUE)), control=list(n=20000, s=50))max=max(bounds[x,]))))

    .geigerwarn(immediate.=TRUE)
    
    ct=list(n=1e4, s=100, w=1, sample.priors=FALSE)
    ct[names(control)]=control
	par=argn(lik)
            
    ## PRIOR SPECIFICATIONS
    # uniform
    unipr=function(a,b) function(x) dunif(x, min=a, max=b, log=TRUE)
    
    # truncated normal (from msm:::dtnorm)
    tnormpr=function(a,b,m,s) {
        denom=pnorm(b, mean=m, sd=s) - pnorm(a, mean=m, sd=s)
        function(x) {
            ret=numeric(length(x))
            if(any(x>=a & x<=b)->ww) {
                ret[ww]=dnorm(x[ww], mean=m, sd=s, log=TRUE)-log(denom)
            }
            if(any(!ww)) ret[!ww]=-Inf
            
            ret
        }
    }
        
    ## PRIORS
    # default
    priors=list(
        z0=unipr(-1e6, 1e6),
        sigsq=tnormpr(-500, 100, m=-500, s=100),
        alpha=tnormpr(-500, 5, m=-500, s=100),
        a=unipr(-10, 10),
        drift=unipr(-100, 100),
        slope=unipr(-100, 100),
        lambda=unipr(0, 1),
        kappa=unipr(0, 1),
        delta=unipr(0, 2.999999),
        SE=tnormpr(-500, 100, m=-500, s=100),
        trns=tnormpr(-500, 100, m=-500, s=100)
    )
    
    # parameter space
    epar=c("z0", "sigsq", "alpha", "a", "drift", "slope", "lambda", "kappa", "delta", "SE", "trns")
    names(epar)=c("nat", "exp", "exp", "nat", "nat", "nat", "nat", "nat", "nat", "exp", "exp")

    for(i in names(priors)){
        xx=match(i, epar)
        attr(priors[[i]], "type")=names(epar)[xx]
    }
    
    priors[names(prior)]=prior
    
    if(!all(par%in%names(prior))) {
        warning(paste("Default prior(s) used for:\n\t", paste(par[!par%in%names(prior)], collapse="\n\t"), sep=""))
    }
    
    if(!all(par%in%names(priors))) {
        flag=paste("Missing prior(s) for:\n\t", paste(missingpar<-par[!par%in%names(priors)], collapse="\n\t"), sep="")
        attb=attributes(lik)
        if(is.null(trns<-attb$trns)) stop(flag)
        if(any(trns)) pp=lapply(par[which(trns)], function(x) priors$trns) else stop(flag)
        names(pp)=par[which(trns)]
        priors=c(priors, pp)
    }
    
    priors=priors[par]
    chk=sapply(priors, function(x) {
        if(!is.null(attributes(x)$type)){
            attributes(x)$type%in%c("exp","nat")
        } else {
            return(FALSE)
        }
    })
    if(any(!chk)) stop("Ensure that priors are given as a list and that each prior has a 'type' attribute that is either 'exp' or 'nat'")

    ## BOUNDS
    bounds=lapply(priors, function(x) {
        ee=environment(x)
        if(is.null(ee$a)) a=-Inf else a=ee$a
        if(is.null(ee$b)) b=Inf else b=ee$b
        return(list(min=a, max=b))
    })
    
    ## REPARAMETERIZED LIK FUNCTION
    typ=sapply(priors, function(x) attributes(x)$type)
    flik=function(p){
        pp=ifelse(typ=="exp", exp(p), p)
        lik(pp)
    }

    getstart=function(bounds){
        fx=function(bound){
            bb=unlist(bound)
            if(any(is.infinite(bb))){
                bb[which(bb==min(bb))]=max(c(-100, min(bb)))
                bb[which(bb==max(bb))]=min(c(100, max(bb)))
            }
            rr=diff(bb)
            as.numeric(min(bb)+rr*runif(1, 0, 1/3)) # get value from lower third of range
        }
        sapply(as.list(bounds), fx)
    }
	
    ## STARTING VALUES
    count=0
    while(1){
        fstart=FALSE
        flag="'start' should be a named vector of starting values with names in argn(lik), returning a finite likelihood"
        
        if(is.null(start)) start=getstart(bounds[par]) else  start=try(ifelse(typ=="exp", log(start), start), silent=TRUE)
        if(inherits(start, "try-error")) stop(flag)
        if(any(is.na(start))) stop(flag)

        lnLc=try(flik(start), silent=TRUE)
        
        if(inherits(lnLc, "try-error")) count=count+1
        if(!is.numeric(lnLc) | is.na(lnLc)) count=count+1 else break()
        if(count==100) stop(flag)
        next()
    }
    
    ## PROPOSAL FREQUENCIES
	if(is.null(proposal)) {
        proposal=structure(rep(1,length(par)), names=par)
    } else {
        if(!setequal(names(proposal), par)){
            stop("'proposal' must be supplied as a named vector of proposal frequencies with names in argn(lik)")
        } else {
            proposal=proposal[par]
        }
    }
	
	tmp=cumsum(proposal)		
	proposal=structure(tmp/max(tmp), names=par)
	
    ## DATA STRUCTURES
	cur<-structure(start, names=par)
	n_prop<-a_prop<-structure(numeric(length(par)), names=par)
	
    cur<-new<-structure(numeric(length(par)), names=par)
	
	kp=seq(0, ct$n, by=ct$s)
	kp[1]=1
	tab=matrix(NA, nrow=length(kp), ncol=length(par)+1)
	rownames(tab)=kp
	colnames(tab)=c(par,"lnL")
    
    tf=ceiling(seq(1,ct$n, length.out=30))

	for(i in 1:ct$n){
		while(1){
			new<-cur
            if(any(is.na(cur))) stop()
			choice=names(proposal)[min(which(runif(1)<proposal))]
			cpar=cur[[choice]]
			if(runif(1)<0.75) {
				cprop=.proposal.multiplier(cpar, ct$w, bounds[[choice]])
			} else {
				cprop=.proposal.slidingwindow(cpar, ct$w, bounds[[choice]])
			}
			ppar=cprop$v
			new[[choice]]=ppar
			lnLn=try(suppressWarnings(flik(new)), silent=TRUE)
            
		
			if(is.finite(lnLn)){
				lnh=cprop$lnHastingsRatio
				lnp=.dlnratio(cpar, ppar, priors[[choice]])
				if(is.finite(lnp) & is.finite(lnh)) break()
			} else {
                if(runif(1)){
                    lnh=0
                    lnp=0
                    lnLn=lnLc
                    new=cur
                }
            }
		}
		r=.proc.lnR(i, choice, lnLc, lnLn, lnp, lnh, heat=1, control=ct)$r
		n_prop[[choice]]<-n_prop[[choice]]+1
		if(runif(1)<r){
			a_prop[[choice]]<-a_prop[[choice]]+1
			cur<-new
			lnLc<-lnLn
		}
        
        if(i%in%tf) {
            if(i==1) cat("|",rep(" ",9),toupper("generations complete"),rep(" ",9),"|","\n")
            xx=sum(tf==i)
            for(j in 1:xx) cat(". ")
            if(i==ct$n) cat("\n")
        }
        
		if(i%in%kp) tab[which(kp==i),]=c(ifelse(typ=="exp",exp(cur),cur), lnLc)
        if(i==max(kp)){
            cat("\n\n",rep(" ",10),toupper(" sampling summary"),"\n")
            df=structure(data.frame(proposed=n_prop, accepted=a_prop, adoptrate=a_prop/n_prop), rownames=par)
            if(any(is.na(df))) df[is.na(df)]=0
            .print.table(df, digits=c(0,0,4), buffer=6)
        }
	}
	return(tab)
}

.TESTING.root.phylo=function(phy, outgroup, taxonomy=NULL){
    ## GENERAL FUNCTION FOR ROOTING (based on outgroup)
    # taxonomy: classification data.frame with 'species' minimally as a rownames
    
	if(!is.null(sys)) {
		sys=cbind(rn=rownames(sys), taxonomy)
		rows=unique(unlist(sapply(outgroup, function(o) which(sys==o, arr.ind=TRUE)[,1])))
		outgroup=rownames(sys)[rows]
		outgroup=outgroup[outgroup%in%phy$tip.label]
	} else {
		if(!all(outgroup%in%phy$tip.label)) stop("Some 'outgroup' appear missing from 'phy'.")
	}
	
	tips=match(outgroup, phy$tip.label)
	node=getMRCA(phy,tips)
	if(node==Ntip(phy)+1){
		node=getMRCA(phy, (1:Ntip(phy))[-tips])
	}
	rooted=root(phy, node=node, resolve.root=TRUE)
	rooted
}





.TESTING.ultrametricize.phylo=function(phy, trim=c("min","max","mean","depth"), depth=NULL){
	
	phy <- reorder(phy)
    n <- length(phy$tip.label)
    n.node <- phy$Nnode
    xx <- numeric(n + n.node)
    for (i in 1:nrow(phy$edge)) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	
	paths=xx[1:n]
	trim=switch(match.arg(trim),
    min = min(paths),
    max = max(paths),
    mean = mean(paths),
    depth = NULL)
	
	if(is.null(trim)) {
		if(!is.null(depth)) trim=depth else stop("'depth' must be supplied if 'trim=depth'")
	}
	
	tol=diff(range(paths))
	
	cat(paste("Detected maximum difference in root-to-tip path lengths of ",tol,"\n",sep=""))
	rsc=function(phy, curdepth, depth) {phy$edge.length=phy$edge.length*(depth/curdepth); phy}
    
	ww=which(phy$edge[,2]<=n)
	phy$edge.length[ww]=phy$edge.length[ww]+(trim-paths[phy$edge[ww,2]])
	if(trim!=depth && !is.null(depth)) {
		phy=rsc(phy, trim, depth)
	}
	
	if(any(phy$edge.length<0)) warning("ultrametricized 'phy' has negative branch lengths")
    
	return(phy)
}

.TESTING_bind.phylo=function(phy, taxonomy){
    ## phy: a 'rank' level phylogeny (tips of 'phy' should be matchable to taxonomy[,rank])
    ## taxonomy: a mapping from genus, family, order (columns in that order); rownames are tips to be added to constraint tree
    ##		-- 'taxonomy' MUST absolutely be ordered from higher exclusivity to lower (e.g., genus to order)
    ## rank: rank at which groups are assumed to be monophyletic (currently for 'family' only)
    ## returns a nodelabeled constraint tree based on 'phy' and 'rank'-level constraints
	
	
	
    #	oo=order(apply(taxonomy, 2, function(x) length(unique(x))),decreasing=TRUE)
    #	if(!all(oo==c(1:ncol(taxonomy)))){
    #		warning("Assuming 'taxonomy' is not from most to least exclusive")
    #		taxonomy=taxonomy[,ncol(taxonomy):1]
    #	}
	taxonomy=as.data.frame(taxonomy, stringsAsFactors=FALSE)
	rank=colnames(taxonomy)[unique(which(taxonomy==phy$tip.label, arr.ind=TRUE)[,"col"])]
	if(length(rank)!=1) stop("tips in 'phy' must occur in a single column of 'taxonomy'")
	
	tax=taxonomy
	ridx=which(colnames(tax)==rank)
	tax=tax[,ridx:ncol(tax)]
	tips=rownames(tax)
	
	original_taxonomy=taxonomy
	
	
	# PRUNE 'rank'-level tree if some taxa unmatched in 'tips'
	if(any(zz<-!phy$tip.label%in%tax[,rank])){
		warning(paste("taxa not represented in 'tips':\n\t", paste(phy$tip.label[zz], collapse="\n\t"), sep=""))
		phy=.drop.tip(phy, phy$tip.label[zz])
	}
	
	tmp=unique(c(phy$tip.label, phy$node.label))
	all_labels=tmp[!is.na(tmp)&tmp!=""]
	exclude=apply(tax, 1, function(x) !any(x%in%all_labels))
	missing_tips=rownames(tax)[exclude]
	if(length(missing_tips)){
		warning(paste("tips missing data in 'taxonomy' and excluded from constraint tree:\n\t", paste(missing_tips, collapse="\n\t"), sep=""))
	}
	tax=tax[!exclude,]
	
	# find tips that have data but whose data not at 'rank' level (plug in deeper in tree)
	at_rank=tax[,rank]%in%phy$tip.label
	deeper_tips=tax[!at_rank,]
	if(nrow(deeper_tips)){
		ww=apply(deeper_tips, 1, function(x) x[[min(which(x%in%all_labels))]])
		for(i in 1:length(ww)){
			nn=.nodefind.phylo(phy, ww[i])
			if(!is.null(nn)){
				tmp=.polytomy.phylo(names(ww[i]))
				phy=bind.tree(phy, tmp, where=nn)
			}
		}
		phy=compute.brlen(phy, method="Grafen")
		warning(paste("tips missing data at 'rank' level in 'taxonomy' but included in constraint tree:\n\t", paste(rownames(deeper_tips), collapse="\n\t"), sep=""))
	}
	phy$node.label=NULL
	tax=as.matrix(original_taxonomy[at_rank,1:ridx])
    if(ridx==1) {
        rownames(tax)=rownames(original_taxonomy)[at_rank]
        colnames(tax)=colnames(original_taxonomy)[ridx]
    }
    tax=as.data.frame(tax)
	if(any(is.na(tax[,rank]))) stop("Corrupted data encountered when checking taxonomy[,rank]")
	
	rsc=function(phy, age=1) {
		ee=phy$edge.length
		ag=max(heights(phy))
		phy$edge.length=ee*(age/ag)
		phy
	}
	## CREATE 'rank'-level subtrees
	mm=min(phy$edge.length)
	ss=split(tax, tax[,rank])
	subtrees=lapply(1:length(ss), function(idx) {
        curnm=names(ss)[idx]
        x=ss[[idx]]
        rnm=rownames(x)
        y=as.matrix(x)
        d=apply(y, 2, function(z) if(all(z=="") | length(unique(z))==1) return(TRUE) else return(FALSE))
        if(any(d)) {
            nm=colnames(y)
            y=as.matrix(y[,-which(d)])
            colnames(y)=nm[-which(d)]
        }
        if(!length(y)){
            cur=.polytomy.phylo(rnm, mm/2)
        } else {
            cur=phylo.lookup(cbind(y,names(ss)[idx]))
            cur=compute.brlen(cur)
        }
        cur$root.edge=0
        cur=rsc(cur, mm/2)
        return(cur)
        
	})
	names(subtrees)=names(ss)
	
	## PASTE in 'rank'-level subtrees
	contree=glomogram.phylo(phy, subtrees)
	contree=compute.brlen(contree)
	contree
}



