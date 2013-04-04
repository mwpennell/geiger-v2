
medusa <- function(phy, richness=NULL, criterion=c("aicc", "aic"), partitions=NA, model=c("mixed", "bd", "yule"), cut=c("both","stem","node"), init=c(r=0.05, epsilon=0.5), ...){
  
	## CHECK ARGUMENTS
    verbose=FALSE
    initialE=init[["epsilon"]]
    initialR=init[["r"]]
    
    mc=.check.parallel()
    num.cores=getOption("cores")
    
    if(!is.na(partitions)){
        flag="'partitions' should either be NA or an integer, specifying the maximum number of piecewise models to consider"
        if(!is.numeric(partitions)) stop(flag)
        if(partitions<=0) stop(flag)
        stop="partitions"
    } else {
        criterion=match.arg(criterion, choices=c("aicc", "aic"), several.ok=FALSE)
        stop="threshold"
    }
    #	stop=match.arg(stop, choices=c("threshold","partitions"), several.ok=FALSE)
    if(stop=="threshold"){
        #        criterion=match.arg(threshold, choices=c("aicc", "aic"), several.ok=FALSE)
        npartitions=0
    } else {
        npartitions=partitions
        criterion="aicc"
    }
	model=match.arg(model, choices=c("mixed", "bd", "yule"), several.ok=FALSE)
	shiftCut=match.arg(cut, choices=c("both", "stem", "node"), several.ok=FALSE)
    
	## PREPARE DATA 
	## Before determining model.limit, prune tree as necessary (from 'taxon' information in 'richness')
    #	if (nexus) phy <- read.nexus(phy)
	if(!any(c("phylo","multiPhylo")%in%class(phy))) stop("'phy' must either be a phylo or multiPhylo object")
    richness=data.frame(richness, stringsAsFactors=FALSE)
	phyData <- .treedata.medusa(phy=phy, richness=richness, ...) ## modified prune.tree.merge.data for multiple trees (jme)
	## END -- jme
	
	# internal function for running a single tree and richness data.frame (jme)
	medusa_runner=function(phy, richness){
		## Limit on number of piecewise models fitted; based on tree size, aicc correction factor,
		## and flavour of model fitted (i.e. # parameters estimated; birth-death or pure-birth)
		model.limit <- .get.max.model.limit(richness=richness, stop=stop, npartitions=npartitions, model=model, verbose=verbose);
		## Determine correct AICc threshold from tree size (based on simulations)
		## Should be used for interpreting model-fit
        N <- Ntip(phy)
		threshold_N=ifelse(stop=="threshold", .threshold.medusa(N), 0)
		#cat("Appropriate AICc threshold for tree of ", N, " tips is: ", threshold, ".\n\n", sep="");
		
		## Store pertinent information: branch times, richness, descendants
		#cat("Preparing data for analysis... ");
		obj <- .make.cache.medusa(phy=phy, richness=richness, mc=mc, num.cores=num.cores);
		#cat("done.\n");
		
		## Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
		pend.nodes <- seq_len(length(phy$tip.label));   # Calculate pendant splits just once, keep track through various models
		int.nodes <- (length(phy$tip.label)+2):max(phy$edge); # Omit root node
		root.node <- length(phy$tip.label) + 1;
		all.nodes <- c(pend.nodes, root.node, int.nodes);
		
		desc <- list(desc.stem=obj$desc.stem, desc.node=obj$desc.node)
		
		z <- z.orig <- obj$z;
		
		## Pre-fit pendant edges so these values need not be re(re(re))calculated; amounts to ~25% of all calculations
		## Will show particular performance gain for edges with many fossil observations
        #		cat("Optimizing parameters for pendant edges... ");
		tips <- NULL;
		# Will always be shiftCut="stem"; if mixed model, keep only best fit and throw out other in .prefit.medusa
		if (mc)
		{
			tips <- mclapply(pend.nodes, .prefit.medusa, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="stem", criterion=criterion, mc.cores=num.cores);
		} else {
			tips <- lapply(pend.nodes, .prefit.medusa, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="stem", criterion=criterion);
		}
        #		cat("done.\n");
		
		## Pre-fit virgin internal nodes; should deliver performance gain for early models, and especially for large trees
		## Remain useful until a spilt is accepted within the clade
		## Need to incorporate cutAtStem here
		#cat("Pre-calculating parameters for virgin internal nodes... ");
		virgin.stem <- list(); virgin.node <- list();
		if (mc)
		{
			if (shiftCut == "stem" || shiftCut == "both")
			{
				virgin.stem <- mclapply(int.nodes, .prefit.medusa, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="stem", criterion=criterion, mc.cores=num.cores);
			}
			if (shiftCut == "node" || shiftCut == "both")
			{
				virgin.node <- mclapply(int.nodes, .prefit.medusa, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="node", criterion=criterion, mc.cores=num.cores);
			}
		} else {
			if (shiftCut == "stem" || shiftCut == "both")
			{
				virgin.stem <- lapply(int.nodes, .prefit.medusa, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="stem", criterion=criterion);
			}
			if (shiftCut == "node" || shiftCut == "both")
			{
				virgin.node <- lapply(int.nodes, .prefit.medusa, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="node", criterion=criterion);
			}
		}
		virgin.nodes <- list(stem=virgin.stem, node=virgin.node);
		#cat("done.\n\n");
		
		prefit <- list(tips=tips, virgin.nodes=virgin.nodes);
		
		## Needed downstream; do not recalculate
		## Gives the number of tips associated with an internal node; determines whether a node is 'virgin' or not
		num.tips <- list()
		if (mc)
		{
			num.tips <- mclapply(all.nodes, function(x) length(obj$tips[[x]]))
		} else {
			num.tips <- lapply(all.nodes, function(x) length(obj$tips[[x]]))
		}
		
		## Fit the base model
		## 'fit' holds current results; useful for initializing subsequent models
		fit <- list();
		if (model == "mixed")
		{
			fit.bd <- .initial.medusa(z=z, initialR=initialR, initialE=initialE, model="bd");
			fit.yule <- .initial.medusa(z=z, initialR=initialR, initialE=initialE, model="yule");
			if (fit.bd[[criterion]] < fit.yule[[criterion]]) {
				fit <- fit.bd;
				fit$model <- "bd";
			} else {
				fit <- fit.yule;
				fit$model <- "yule";
			}
		} else if (model == "bd") {
			fit <- .initial.medusa(z=z, initialR=initialR, initialE=initialE, model="bd");
			fit$model <- "bd";
		} else {
			fit <- .initial.medusa(z=z, initialR=initialR, initialE=initialE, model="yule");
			fit$model <- "yule";
		}
        
		models <- list(fit);
		zz <- list(z)
        
		if (stop == "partitions")
		{

            #			cat("Step 1 (of ", model.limit, "): best likelihood = ", models[[1]]$lnLik, "; AICc = ", models[[1]]$aicc, "; model = ", models[[1]]$model, "\n", sep="");
			for (i in seq_len(model.limit-1))
			{
				node.list <- all.nodes[-fit$split.at];
				if (mc)  # parallel (i.e. multithreaded) processing. No GUI, and not at all on Windows
				{
					res <- mclapply(node.list, .update.fit.medusa, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, shiftCut=shiftCut, mc.cores=num.cores);
				} else {
					res <- lapply(node.list, .update.fit.medusa, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, shiftCut=shiftCut);
				}
				# Select model with best score according to the specific criterion employed (default aicc)
				best <- which.min(unlist(lapply(res, "[[", criterion)));
				models <- c(models, res[best]);
				fit <- res[[best]];   # keep track of '$split.at' i.e. nodes already considered
				
				z <- zz[[i+1]] <- .split.z.at.node.medusa(node=node.list[best], z=z, desc=desc, shiftCut=fit$cut.at)$z;
				
                #				cat("Step ", i+1, " (of ", model.limit, "): best likelihood = ", round(models[[i+1]]$lnLik, digits=7), "; AICc = ", models[[i+1]]$aicc, "; break at node ", models[[i+1]]$split.at[i+1], "; model=", models[[i+1]]$model, "; cut=", models[[i+1]]$cut.at, "\n", sep="");
			}
		} else if (stop == "threshold") {
			i <- 1;
			done <- FALSE;
            #		cat("Step 1: best likelihood = ", models[[1]]$lnLik, "; AICc = ", models[[1]]$aicc, "; model = ", models[[1]]$model, "\n", sep="");
			while (!done & i < model.limit)
			{
				node.list <- all.nodes[-fit$split.at];
				if (mc)  # parallel (i.e. multithreaded) processing. No GUI, and not at all on Windows
				{
					res <- mclapply(node.list, .update.fit.medusa, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, shiftCut=shiftCut, mc.cores=num.cores);
				} else {
					res <- lapply(node.list, .update.fit.medusa, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, shiftCut=shiftCut);
				}
				# Select model with best score according to the specific criterion employed (default aicc)
				best <- which.min(unlist(lapply(res, "[[", criterion)));
				# Compare last accepted model to current best model
				if (as.numeric(models[[length(models)]][criterion]) - as.numeric(res[[best]][criterion]) < threshold_N)
				{
                    #				cat("\nNo significant increase in ", criterion, " score. Disregarding subsequent piecewise models.\n", sep="");
					done <- TRUE;
					break;
				}
				models <- c(models, res[best]);
				fit <- res[[best]];   # keep track of '$split.at' i.e. nodes already considered
				
				z <- zz[[i+1]] <- .split.z.at.node.medusa(node=node.list[best], z=z, desc=desc, shiftCut=fit$cut.at)$z;
				
                #				cat("Step ", i+1, ": best likelihood = ", models[[i+1]]$lnLik, "; AICc = ", models[[i+1]]$aicc, "; break at node ", models[[i+1]]$split.at[i+1], "; model=", models[[i+1]]$model, "; cut=", models[[i+1]]$cut.at, "\n", sep="");
				i <- i+1;
			}
		}
		
		modelSummary <- .summary.modelfit.medusa(models=models, phy=phy, threshold=threshold_N)
        
		if (verbose)
		{
			cat("\n", "Model fit summary:", "\n\n", sep="");
			print(modelSummary);
			if (threshold_N > 0)
			{
				cat("\nAIC weights are not reported, as they are meaningless when using a threshold criterion.\n")
			}
		}
		
        
		## SUMMARY
		zSummary=function(id){
            model.id=id
            z=zz[[model.id]]
			opt.model = data.frame(cut=modelSummary$cut[1:model.id], split=modelSummary$split[1:model.id], models[[model.id]]$par, lnLp=models[[model.id]]$lnLik.part, stringsAsFactors=FALSE)
			break.pts = opt.model$split
			cut.at = opt.model$cut
			rownames(z)=phy$hash[z[,"dec"]] # identify identical edges among trees
			
			# collect shift times
			internal=is.na(z[,"n.t"])
			res=numeric(nrow(z))
			res[]=NA
			check=!is.na(break.pts)
			cut.at=cut.at[check]
			break.pts=break.pts[check]
			if(length(break.pts)){
				for(i in 1:length(break.pts)){
					idx=which(z[,"dec"]==break.pts[i])
					if(cut.at[i]=="node" && internal[idx]) {
						res[idx]=z[idx,"t.1"]
					} else if(cut.at[i]=="stem" || !internal[idx]){
						res[idx]=z[idx,"t.0"]
					} else {
						stop("'cut' should either be stem or node")
					}
				}
			}
			anc=match(z[,"anc"],z[,"dec"])
			z=cbind(z, shift=match(z[,"dec"], opt.model$split), t.shift=res)
			z=cbind(z,r=opt.model[z[,"partition"],"r"],epsilon=opt.model[z[,"partition"],"epsilon"])
			z=cbind(z, ancestral.r=z[anc,"r"], ancestral.epsilon=z[anc,"epsilon"])
			return(z)
		}
        
        model.id=nrow(modelSummary)
		z.summary<-zSummary(id=model.id)
		
        control=list(stop=stop, threshold=structure(threshold_N, names=criterion), partitions=npartitions)
		results <- list(control=control, cache=list(desc=desc, phy=phy, richness=richness), models=models, summary=modelSummary, FUN=zSummary);
		return(results)
	}
	
	if("multiPhylo"%in%class(phy)) { ## deal with multiple trees
		res=lapply(phyData, function(x) medusa_runner(x$phy, x$richness))
		names(res)=names(phy)
		class(res)=c("multimedusaRAW",class(res))
	} else {
		phyData$phy=hashes.phylo(phyData$phy)
		res=medusa_runner(phyData$phy, phyData$richness)
		class(res)=c("medusaRAW",class(res))
	}
	
	return(res)
}