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
        SE=tnormpr(-500, 100, m=-500, s=100)
    )
    
    # parameter space
    epar=c("z0", "sigsq", "alpha", "a", "drift", "slope", "lambda", "kappa", "delta", "SE")
    names(epar)=c("nat", "exp", "exp", "nat", "nat", "nat", "nat", "nat", "nat", "exp")

    for(i in names(priors)){
        xx=match(i, epar)
        attr(priors[[i]], "type")=names(epar)[xx]
    }
    
    priors[names(prior)]=prior
    
    if(!all(par%in%names(prior))) {
        warning(paste("Default prior(s) used for:\n\t", paste(par[!par%in%names(prior)], collapse="\n\t"), sep=""))
    }
    
    if(!all(par%in%names(priors))) {
        stop(paste("Missing prior(s) for:\n\t", paste(par[!par%in%names(priors)], collapse="\n\t"), sep=""))
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

    
    ## REPARAMETERIZED LIK FUNCTION
    typ=sapply(priors, function(x) attributes(x)$type)
    flik=function(p){
        pp=ifelse(typ=="exp", exp(p), p)
        lik(pp)
    }

	
    ## STARTING VALUES
    fstart=FALSE
    if(is.null(start)) {
        start=structure(rep(1, length(par)), names=par)
        fstart=TRUE
    }
    if(is.null(names(start)) & length(start)!=length(par)) stop("'start' must be supplied as a named vector of starting values with names in argn(lik)")
    if(!is.null(names(start)) & all(names(start)%in%par)) start=start[par] else if(!length(par)==length(start)) stop("'start' appears misspecified")
    
    start=ifelse(typ=="exp", log(start), start)
    if(!is.numeric(flik(start))) {
        if(fstart) warning("try supplying 'start' as a named vector of starting values for each parameter in argn(lik)")
        stop("'start' must return a numeric likelihood for 'lik'")
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
	lnLc=flik(start)
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
				cprop=.proposal.multiplier(cpar, ct$w, list(min=-Inf, max=Inf))
			} else {
				cprop=.proposal.slidingwindow(cpar, ct$w, list(min=-Inf, max=Inf))
			}
			ppar=cprop$v
			new[[choice]]=ppar
			lnLn=flik(new)
		
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
