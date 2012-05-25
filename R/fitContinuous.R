`fitContinuous` <-
function(phy, data, data.names=NULL, model=c("BM", "OU", "lambda", "kappa", "delta", "EB", "white", "trend"), bounds=NULL,  meserr=NULL)
{
	
	# sort is T because sub-functions assume data are in
	# this particular order
	
	model<-match.arg(model)
	
	td<-treedata(phy, data, data.names, sort=T)

	ntax=length(td$phy$tip.label)

	if(is.null(meserr)) {
		me=td$data
		me[]=0
		meserr=me	
	} else if(length(meserr)==1) {
		me=td$data
		me[]=meserr
		meserr=me
	} else if(is.vector(meserr)) {
		if(!is.null(names(meserr))) {
			o<-match(rownames(td$data), names(meserr))
			if(length(o)!=ntax) stop("meserr is missing some taxa from the tree")
			meserr<-as.matrix(meserr[o,])
		} else {
			if(length(meserr)!=ntax) stop("No taxon names in meserr, and the number of taxa does not match the tree")
			me<-td$data
			me[]=meserr
			meserr=me
		}
	} else {
		if(!is.null(rownames(meserr))) {
			o<-match(rownames(td$data), rownames(meserr))
			meserr=meserr[o,]
		} else {
			if(sum(dim(meserr)!=dim(td$data))!=0)
				stop("No taxon names in meserr, and the number of taxa does not match the tree")
			print("No names in meserr; assuming that taxa are in the same order as tree")	
		}
	}

	#--------------------------------
    #---    PREPARE DATA LIST     ---
    #--------------------------------
	ds			<- list()
   		ds$tree 		<- td$phy          # TIP data 
    #--------------------------------
    #--- SET MODEL SPECIFICATIONS ---
    #--------------------------------
    cat("Fitting ", model, "model:\n")
    #-----------------------------
    #---  SET PARAMETER BOUNDS ---
    #-----------------------------
    #---- DEFAULT BOUNDS
    bounds.default			 <- matrix(c(0.00000001, 20, 0.0000001,1, 0.000001, 1, 0.00001, 2.999999, 0.0000001, 50, -3, 0, 0.0000000001, 100, -100, 100), nrow=8, ncol=2, byrow=TRUE)
    rownames(bounds.default) <- c("beta", "lambda", "kappa", "delta", "alpha", "a", "nv", "mu");
    colnames(bounds.default) <- c("min", "max")

 	#---- USER DEFINED PARAMETER BOUNDS
 	if (is.null(bounds)) {
 		bounds <- bounds.default       # USE DEFAULTS
 	}else{
 		if (class(bounds)!="list"){
 			stop("Please specify user defined parameter bounds as a list()")
 		}else{
 			specified   <- !c(is.null(bounds$beta), is.null(bounds$lambda), 
 							  is.null(bounds$kappa), is.null(bounds$delta),  is.null(bounds$alpha), is.null(bounds$a),
 							  is.null(bounds$nv), is.null(bounds$mu)
 							  )
 			bounds.user <- matrix(c(bounds$beta, bounds$lambda, bounds$kappa, bounds$delta, bounds$alpha, bounds$a, bounds$nv, bounds$mu), 
 								  nrow=sum(specified), ncol=2, byrow=TRUE
 								  )
 			rownames(bounds.user) <- c("beta", "lambda", "kappa", "delta", "alpha", "a",  "nv", "mu")[specified]
   	 		colnames(bounds.user) <- c("min", "max")
  
   	 		#----  SET FINAL SEARCH BOUNDS
 			bounds <- bounds.default
 			bounds[specified,] <- bounds.user     # Final Bounds
   		} # END if list
   	}  # END user bound if loop
   	#--------------------------------
    #---   APPEND MODEL SETTINGS  ---
    #--------------------------------
  	ds$bounds <- data.frame(t(bounds))
	print(ds$bounds)
  	ds$model  <- model
  	#--------------------------------
    #---        FIT MODEL         ---
    #--------------------------------
    result<-list()
    for(i in 1:ncol(td$data)) {
    	ds$data=td$data[,i]
    	ds$meserr=meserr[,i]
  		result[[i]]<-fitContinuousModel(ds, print=print)
  		if(!is.null(colnames(td$data))) names(result)[i]<-colnames(td$data)[i] else names(result)[i]<-paste("Trait", i, sep="")

  	}
  	result
}


`fitContinuousModel` <-
function(ds, print=TRUE)
{
	bounds 	<- ds$bounds
	model 	<- ds$model
	n 		<- length(ds$data)
	print(bounds)
	#----- MINIMIZE NEGATIVE LOG LIKELIHOOD
	
	beta.start<-var(ds$data)/max(branching.times(ds$tree))


	out         <- NULL
	
	y			<- ds$data				# TIP data
	tree		<- ds$tree			# Tree
	meserr		<- ds$meserr
	n			<- length(y)
	

	#----------------------------------
	#-----       DEFAULT FIT      -----
	#----------------------------------
	if (model=="BM") {
		k<-2
	
		vcv<-vcv.phylo(tree)

		start=log(beta.start)
		lower=log(bounds[1,"beta"])
		upper=log(bounds[2,"beta"])
		
		foo<-function(x) {
			vv<-exp(x)*vcv
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			-dmvnorm(y, mu, vv, log=T)
		}
		
		o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		sigSquared<-exp(o$par)
		vv<-exp(sigSquared)*vcv
		diag(vv)<-diag(vv)+meserr^2
		zBar<-phylogMean(vv, y)[1,1]
		results<-list(lnl=-o$value, zBar=zBar, sigSquared= sigSquared)

	#----------------------------------
	#-----       LAMBDA ONLY      -----
	#----------------------------------
	} else if (model=="lambda"){
		k<-3
		
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","lambda")])
		upper=log(bounds[2,c("beta","lambda")])
		
		
		foo<-function(x) {


			vcv<-vcv.phylo(tree)

			index			<-	matrix(TRUE, n,n)
			diag(index)		<- FALSE
			vcv[index] 	<- vcv[index]*exp(x[2])
			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
			
		results<-list(lnl=-o$value, beta= exp(o$par[1]), lambda=exp(o$par[2]))

	#----------------------------------
	#-----        KAPPA ONLY      -----
	#----------------------------------
	} else if (model=="kappa"){
		k<-3
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","kappa")])
		upper=log(bounds[2,c("beta","kappa")])
				
		
		foo<-function(x) {

			t<-kappaTree(tree, kappa=exp(x[2]))
			vcv<-vcv.phylo(t)

			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		results<-list(lnl=-o$value, beta= exp(o$par[1]), kappa=exp(o$par[2]))


	#----------------------------------
	#-----        DELTA ONLY      -----
	#----------------------------------	
	} else if (model=="delta"){
		
		k<-3
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","delta")])
		upper=log(bounds[2,c("beta","delta")])
		
		foo<-function(x) {

			t<-deltaTree(tree, delta=exp(x[2]))
			vcv<-vcv.phylo(t)

			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		results<-list(lnl=-o$value, beta= exp(o$par[1]), delta=exp(o$par[2]))	
	#----------------------------------
	#-----        WHITE NOISE     -----
	#----------------------------------	
	} else if (model=="white"){
		
		k<-2
		start=c(mean(y), log(var(y)))
		lower=c(-Inf, log(bounds[1,"nv"]))
		upper=c(Inf, log(bounds[2, "nv"]))
		
		lnl.noise<- function (p, x, se)
		# p is the vector of parameters, tree is not needed
		# x and se are trait means and std errors
		{
  			## prep parameters
  			root<- p[1]	# trait value of root of tree (also optimum)
  			vs<- exp(p[2])		# white noise variance 
  			n<- length(x)
  			VV<- diag(vs, nrow=n)	
  			diag(VV)<- diag(VV) + se^2	# add sampling variance to diagonal
    
  			## logl
  			M<- rep(root,n)
  			-dmvnorm(x, M, VV, log=TRUE)
		}
		
		o<- optim(start, fn=lnl.noise, x=y, se=meserr, lower=lower, upper=upper, method="L")
		
		results<-list(lnl=-o$value, mean= o$par[1], nv=exp(o$par[2]))	
	#----------------------------------
	#-----        TREND           -----
	#----------------------------------	
	} else if (model=="trend"){
		
		k<-3
		vcv<- vcv.phylo(tree)
  		ww<- lm(y ~ diag(vcv))
  		p0<- c(phylogMean(vcv, y), var(y)/max(branching.times(tree)), coef(ww)[2])
		if(is.na(p0[3])) {
			p0[3]<-0
			if(is.ultrametric(tree))
				cat("WARNING: Cannot estimate a trend with an ultrametric tree; lnl will be the same as the BM model")
		}
		lower=c(-Inf, log(bounds[1,"beta"]), bounds[1,"mu"])
		upper=c(Inf, log(bounds[2,"beta"]), bounds[2,"mu"])
		

		lnl.BMtrend<- function(p, vcv, x, se)
		# p is vector of parameters, tr is tree
		# x and se are vectors of trait means and standard errors
		{
  			## prep parameters
  			root<- p[1]	# trait value of root of tree
  			vs<- exp(p[2])		# BM variance 
  			ms<- p[3]		# BM trend
  			VV<- vs*vcv	
  			diag(VV)<- diag(VV) + se^2	# add sampling variance to diagonal
  			n<- length(x)
  
 			## logl
  			M<- root+ ms*diag(vcv)
  			- dmvnorm(x, M, VV, log=TRUE)
  		}

		o<- optim(p0, fn=lnl.BMtrend, vcv=vcv, x=y, se=meserr, lower=lower, upper=upper, method="L")
		names(o$par)<-NULL
		results<-list(lnl=-o$value, mean= o$par[1], beta=exp(o$par[2]), mu=o$par[3])		
	#----------------------------------
	#-----        ALPHA ONLY      -----
	#----------------------------------			
	} else if (model=="OU"){
	## modified 12 dec 07 to call ouMatrix(x) instead of vcv.phylo(ouTree(x))

		k<-3
		
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","alpha")])
		upper=log(bounds[2,c("beta","alpha")])
	
		print(lower)
		print(upper)
	
		vcvOrig<-vcv.phylo(tree)
		foo<-function(x) {
			vcv <- ouMatrix(vcvOrig, exp(x[2]))
			
			## t<-ouTree(tree, exp(x[2]))
			##vcv<-vcv.phylo(t)
			
			vv<-exp(x[1])*vcv
			diag(vv)<-diag(vv)+meserr^2
			
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		
		outTries<-list()
		
		# First one: try near BM solution
		start=c(log(beta.start), -50)
		outTries[[1]]<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		# Second one: try one with very strong constraints
		tv<-var(y)
		start=log(c(tv*2000, 1000))
		outTries[[2]]<-optim(foo, p=start, lower=lower, upper=upper, method="L")
	

	
		# Try ten random ones
		for(i in 1:10){
			while(1) {

				ll=c(runif(2, min=-20, max=-1))
				uu=ll+runif(2, min=0, max=10)
				start=c(runif(1, min=ll[1], max=uu[1]), runif(1, min=ll[2], max=uu[2]))
				te<-try(outTries[[i+2]]<-optim(foo, p=start, lower=lower, upper=upper, method="L"), silent=T)
				if(class(te)!="try-error") break
				}
				
		}
		
		# Try range of alphas
		atry<- -5:4
		stry<- log(tv*2*exp(atry))
		for(i in 1:10){
			while(1) {

				ll=c(-20, -20)
				uu=c(10, 10)
				start=c(stry[i], atry[i])
				te<-try(outTries[[i+12]]<-optim(foo, p=start, lower=lower, upper=upper, method="L"), silent=T)
				if(class(te)!="try-error") break
				}
				
		}
		
		
		
		ntries<-22
		ltry<-numeric(ntries)
		lsol<-matrix(nrow= ntries, ncol=2)
		for(j in 1:ntries) {
				ltry[j]<-outTries[[j]]$value
				lsol[j,]<-exp(outTries[[j]]$par)
			}

		ltd<-ltry-min(ltry)
		b<-min(which(ltry==min(ltry)))

		gc<-which(ltd<0.01)
		us<-lsol[gc,1]
		usc<-sum((us-min(us))>0.01)			
		out<-outTries[[b[1]]]	
		if(usc>1) {out$message="Warning: likelihood surface is flat."}
			
		if(out$convergence!=0) {out$message="Warning: may not have converged to a proper solution."}

		sigSquared<-exp(out$par[1])
		alpha=exp(out$par[2])
		
		vcv <- ouMatrix(vcvOrig, alpha)
		vv<-sigSquared*vcv
		diag(vv)<-diag(vv)+meserr^2
			
		zBar <-phylogMean(vv, y)[1,1]

		
		
		results<-list(lnl=-out$value, zBar=zBar, sigSquared= sigSquared, alpha=alpha, convergence=out$convergence, message=out$message, k=k)


	#----------------------------------
	#-----        EB ONLY      -----
	#----------------------------------	
	} else if(model=="EB"){

		k<-3
		start=c(log(beta.start), 0.01)
		lower=c(log(bounds[1,"beta"]),bounds[1,"a"])
		upper=c(log(bounds[2,"beta"]),bounds[2,"a"])
		
		foo<-function(x) {
			t<-exponentialchangeTree(tree, a=x[2])

			vcv<-vcv.phylo(t)
			
			vv<-exp(x[1])*vcv
			diag(vv)<-diag(vv)+meserr^2
			
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		sigSquared<-exp(o$par[1])
		r<-o$par[2]
		
		t<-exponentialchangeTree(tree, a=r)

		vcv<-vcv.phylo(t)
			
		vv<-sigSquared*vcv
		diag(vv)<-diag(vv)+meserr^2
			
		zBar<-phylogMean(vv, y)[1,1]

				
		
		results<-list(lnl=-o$value, zBar=zBar, sigSquared= sigSquared, r=r)

	}
	results$aic<-2*k-2*results$lnl
	results$aicc<-2*k*(n-1)/(n-k-2)-2*results$lnl
	results$k<-k
	return(results) 
}



phylogMean<-function(phyvcv, data) 
{
	o<-rep(1, length(data))
	ci<-solve(phyvcv)
	
	m1<-solve(t(o) %*% ci %*% o)
	m2<-t(o) %*% ci %*% data
	
	return(m1 %*% m2)
	
	}
	
ouMatrix <- function(vcvMatrix, alpha) 
{
## follows Hansen 1997; does not assume ultrametricity (AH 12 dec 07)
## vectorized by LJH
  vcvDiag<-diag(vcvMatrix)
  diagi<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag))
  diagj<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag), byrow=T)

  Tij = diagi + diagj - (2 * vcvMatrix)
    
  vcvRescaled = (1 / (2 * alpha)) * exp(-alpha * Tij) * (1 - exp(-2 * alpha * vcvMatrix))
  return(vcvRescaled) 
}


    	