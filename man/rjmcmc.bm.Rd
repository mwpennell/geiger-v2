\name{rjmcmc.bm}
\alias{rjmcmc.bm}
\title{Bayesian sampling of shifts in trait evolution: relaxed Brownian motion}
\description{Implements reversible-jump Markov chain Monte Carlo sampling for trait evolutionary models}
\usage{
rjmcmc.bm(phy, dat, SE=NA, ngen = 50000, sample.freq = 100, type = c("jump-rbm", "rbm", "jump-bm", "bm"), ...)
}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{dat}{a named vector of continuous trait values, associated with each species in \code{phy}}
  \item{SE}{a named vector of standard errors for each trait value; applied to all trait values if given a single value}
  \item{ngen}{number of sampling generations}
  \item{sample.freq}{frequency with which Markov samples are retained (e.g., \code{sample.freq=10} retains every tenth sample in the chain)}
  \item{type}{the class of model to use (see \bold{Details})}
  \item{...}{arguments passed to \code{\link[geiger]{cache}}}
}

\details{
Implemented is an MCMC sampler for a general model of Brownian motion, which in the full model (\code{jump-rbm}) allows relaxed local clocks and also a point process of pulses 
in evolutionary rate along individual branches. Restricted models include global-rate Brownian motion (\code{jump-rbm}), relaxed-rates Brownian motion (\code{rbm}), and models including
jumps but a single rate of diffusion across the tree (\code{jump-bm}). 

Where applicable, posterior estimates of \code{shifts} between local rates, estimates of the \code{rates} themselves, and inferred \code{jumps} (or pulses) are provided as output. 
Estimates are stored as an MCMC-generations-by-branches matrix (see \bold{Examples}), and branches are uniquely labeled by a cryptographic function to ensure comparability amongst trees 
differing in topology (see \code{\link[digest]{digest}}). 

Note that default settings (as the user assumes if nothing is specified in \code{...}) provide absolutely no guarantee of the chain achieving convergence.  The user is 
emphatically encouraged to supply informed arguments for what are the most critical aspects of this MCMC sampler (see \code{\link[geiger]{cache}} for more information on 
permissible modifications to the MCMC sampler). Finding reasonable run parameters will likely require much trial and error.  Run diagnosis and inspection of chain mixing is 
facilitated by the R-package \pkg{coda} or by the Java application, Tracer (\url{http://tree.bio.ed.ac.uk/software/tracer/}).

In the \bold{Examples} below, do \emph{not} expect such short chains to reach stationarity!

}
\value{
After a run has completed, acceptance rates for the primary proposal mechanisms are printed to the console, along with 
settings of control parameters for the run (see \code{\link[geiger]{cache}}).  

Posterior results are written to several files within a base directory, the contents of which are as follows:
  \item{\bold{log} }{is a logfile including the following for each Markov chain: the generations at which samples were retained (\code{state}), the \code{min}, \code{max}, and \code{median} rate of the 
  diffusion process across the tree, the number of evolutionary pulses (\code{jumps}) along single branches, the variance associated with the jump process (\code{jumpvar}), the \code{root} state, 
  and the likelihood (\code{lnL}) and prior (\code{lnLp}) of the model at sampled generations.}
  
  \item{\bold{rda} }{is a compressed R object, which stores branchwise estimates of the jump and diffusion processes. In order to be interpretable, the \code{rda} file should be processed by the function 
  \code{\link[geiger]{load.rjmcmc}}.  The package \pkg{coda} can be used within R for run diagnostics on the processed output (see, e.g., \code{\link[coda]{heidel.diag}} and \code{\link[coda]{autocorr}}.}
}

\author{JM Eastman, LJ Harmon, and AL Hipp}

\seealso{
\code{\link[geiger]{load.rjmcmc}}
}

\examples{
## generate data
jumptarget=10

## DEFINE FUNCTIONS
# jump-diffusion simulator
levylite=function(phy, alpha=0, sigmasq.brown=0.01, sigma.jump=0.2, lambda.jump=0.2){
	phy=reorder(phy)
	cedges=cumsum(phy$edge.length)
	tmax=sum(phy$edge.length)
	nn=phy$edge[,2]
	cs=c()
	tt=c()
	jumps=0
	if(lambda.jump>0){
		while(1){
			dt=rexp(1,lambda.jump)
			tt=c(tt,dt)
			cs=cumsum(tt)
			if(cs[length(cs)]>tmax){
				cs=cs[-length(cs)]
				jumps=length(cs)
				break()
			}
		}		
	}

	jump.edges=rep(0, length(nn))
	if(jumps>0) {
		for(j in 1:length(cs)){
			tmp=min(which(cs[j]<cedges))
			jump.edges[tmp]=jump.edges[tmp]+1
		}
	} 
	
	hist=as.data.frame(matrix(cbind(phy$edge, phy$edge.length, NA, jump.edges, NA, NA, NA), ncol=8))
	names(hist)=c("ancestor","descendant","edge","phenotype","jumps","time", "effect_Brown", "effect_Jump")

	root=Ntip(phy)+1
	for(i in 1:nrow(hist)){
		start=ifelse(hist$ancestor[i]==root, alpha, hist$phenotype[which(hist$descendant==hist$ancestor[i])])
		stime=ifelse(hist$ancestor[i]==root, 0, hist$time[which(hist$descendant==hist$ancestor[i])])
		t=hist$edge[i]
		hist$phenotype[i]=start+(Beffect<-rnorm(1, mean=0, sd=sqrt(sigmasq.brown*t)))
		hist$effect_Brown[i]=abs(Beffect^2)/hist$edge[i]
		hist$time[i]=stime+hist$edge[i]
		if((jmp<-hist$jumps[i])>0){
			hist$phenotype[i]=hist$phenotype[i]+(Jeffect<-sum(rnorm(jmp, mean=0, sd=sigma.jump)))
		} else {
			Jeffect=0
		}
		hist$effect_Jump[i]=abs(Jeffect^2)/hist$edge[i]
	}
	
	return(hist)
}

# plotting
traitgram=function(phy, hist, alpha=0, cex.node=2, scl=2, ...){
	mm=max(abs(alpha-unlist(hist$phenotype)))
	root=Ntip(phy)+1
	plot(x=NULL, y=NULL, xlim=range(pretty(c(0,max(hist$time)))), ylim=range(pretty(c(-mm+alpha, mm+alpha))), bty="n", xlab="time", ylab="phenotypic value")
	for(i in 1:nrow(hist)) {
		start=ifelse(hist$ancestor[i]==root, alpha, hist$phenotype[which(hist$descendant==hist$ancestor[i])])
		stime=ifelse(hist$ancestor[i]==root, 0, hist$time[which(hist$descendant==hist$ancestor[i])])
		
		end=hist$phenotype[i]
		etime=hist$time[i]
		lines(c(stime,etime),c(start,end),col=geiger:::.transparency("gray25",0.75),...)		
	}
	points(hist$time,hist$phenotype,bg=geiger:::.transparency("white",0.75),pch=21,cex=ifelse(hist$descendant<=Ntip(phy), cex.node, cex.node/scl))	
	points(0,alpha,bg=geiger:::.transparency("white",0.75),pch=21,cex=cex.node/scl)		
}

# wrapper for levylite
jumper=function(phy, alpha=0, sigmasq.brown=0.001, sigmasq.jump=0.1, jumps=jumptarget){
	lambda=jumps/sum(phy$edge.length)
	hist=levylite(phy, alpha=alpha, sigmasq.brown=sigmasq.brown, sigma.jump=sqrt(sigmasq.jump), lambda.jump=lambda)
	dat=hist$phenotype[hist$descendant<=Ntip(phy)]
	names(dat)=phy$tip.label[hist$descendant[hist$descendant<=Ntip(phy)]]
	edges=hist$descendant[hist$jumps>0]
	return(list(hist=hist, dat=dat, phy=phy))
}

# scale branches of a subtree
scaler=function(phy, scl=64, min=4, ...){
	while(1) {	
		# find an internal edge
		anc=geiger:::.get.desc.of.node(Ntip(phy)+1,phy)
		branches=phy$edge[,2]
		ss=sapply(anc, function(x) length(geiger:::.get.descendants.of.node(x, phy, tips=FALSE)))
		if(all(ss<min)) stop("'min' is too large")
		branches=branches[branches>Ntip(phy) & branches!=anc]
		branch=branches[sample(1:length(branches),1)]
		desc=geiger:::.get.descendants.of.node(branch,phy)
		if(length(desc)>=min) break()
	}
	rphy=phy
	rphy$edge.length[match(desc,phy$edge[,2])]=phy$edge.length[match(desc,phy$edge[,2])]*scl
	e=numeric(nrow(phy$edge))
	e[match(c(branch,desc),phy$edge[,2])]=1
	cols=c("red","gray")
	dev.new()
	plot(phy,edge.col=ifelse(e==1,cols[1],cols[2]), edge.width=2, ...) 
	
	mtext("expected pattern of rates")

	rphy
}



## GENERATE DATA: jump-diffusion
phy=ladderize(rcoal(200), right=FALSE)

defpar=par(no.readonly=TRUE)
tmp=jumper(phy)
dat=tmp$dat
hist=tmp$hist
hist$scl=hist$effect_Jump/hist$effect_Brown
hist$scl=hist$effect_Jump

hist$cex=(hist$scl-min(hist$scl))/(max(hist$scl)-min(hist$scl))

r=paste(sample(letters,9,replace=TRUE),collapse="")

traitgram(phy, hist, alpha=0)

## RUN ANALYSIS
rjmcmc.bm(phy,dat, prop.width=1.5, ngen=20000, sample.freq=500, filebase=r, simple.start=TRUE, type="jump-bm")
outdir=paste("jump-BM", r, sep=".")
ps=load.rjmcmc(outdir)
dev.new()
plot(x=ps, par="shifts", burnin=0.25, legend=TRUE, show.tip=FALSE)

dev.new()
plot(x=ps, par="jumps", burnin=0.25, legend=FALSE, show.tip=FALSE, type="fan", edge.width=2)

mm=match(phy$edge[,2],hist$descendant)
hist=hist[mm,]
edgelabels.auteur(text=NULL, pch=21, cex=4*asin(sqrt(hist$cex)), bg=NA, col=ifelse(hist$cex>0, 1, NA), lty=2)

par(defpar)
dev.new()
autocorr.plot(ps$log, ask=dev.interactive())
plot(ps$log, ask=dev.interactive())


## GENERATE DATA: multi-rate diffusion
scl=scaler(phy, show.tip=FALSE)
dat=rTraitCont(scl)

## RUN ANALYSIS
rjmcmc.bm(phy, dat, prop.width=1.5, ngen=20000, sample.freq=500, filebase=r, simple.start=TRUE, type="rbm")
outdir=paste("relaxedBM", r, sep=".")
ps=load.rjmcmc(outdir)
dev.new()
plot(x=ps, par="shifts", burnin=0.25, legend=TRUE, show.tip=FALSE)

## PASTE UNCOMMENTED FOLLOWING LINE TO DROP DIRECTORIES CREATED BY RJMCMC
 # unlink(dir(pattern=paste(r)),recursive=TRUE)

}
