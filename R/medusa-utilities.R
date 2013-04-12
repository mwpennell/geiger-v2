## Function to prune tree using 'richness' information, assumed to have minimally two columns, "taxon" and "n.taxa"
##   Perhaps relax on these column names, may cause too many problems
## May also include 'exemplar' column; in that case, rename relevant tip.label before pruning.
#prune.tree.merge.data
.treedata.medusa <- function (phy, richness=NULL, ...)
{
	
	## MODIFIED -- jme
	if("multiPhylo"%in%class(phy)) { ## deal with multiple trees
		res=lapply(phy, .treedata.medusa, richness)
		tips=c()
		for(i in 1:length(res)) tips=union(tips, res[[i]]$phy$tip.label)
		nn=sapply(res, function(x) Ntip(x$phy))
		if(!all(nn==length(tips))) stop("'phy' appears to have trees with incompletely overlapping sets of tips")
		trees=lapply(res, "[[", "phy")
		class(trees)="multiPhylo"
		trees=hashes.phylo(phy=trees, tips=tips)  ## FROM PHYLO package [returns trees with $hash object -- a 'key' for each split]
		for(i in 1:length(res)) res[[i]]$phy=trees[[i]]
		return(res)
	}
	
	if(is.null(richness)) richness <- data.frame(taxon=phy$tip.label, n.taxa=1)
    
	## END -- jme
	
    # Rename exemplar taxa with taxon name in richness file
	if (!is.null(richness$exemplar))
	{
        # Change relevant tip.labels in phy; individual 'exemplar' may be NA, use original tip.label.
        # Ordering in richness file should NOT be assumed to match order of tip.labels
		i.na <- is.na(richness$exemplar);
		phy$tip.label[match(richness$exemplar[!i.na], phy$tip.label)] <- as.character(richness$taxon[!i.na]);
	}
	
    # make sure things are in the correct order and of correct format
	if (length(richness[1,]) == 2)
	{
		if (colnames(richness)[1] != "taxon" || colnames(richness)[2] != "n.taxa")
		{
			if (class(richness[,1]) == "factor" & class(richness[,2]) == "integer")
			{
				colnames(richness) = c("taxon", "n.taxa");
			} else if (class(richness[,1]) == "integer" & class(richness[,2]) == "factor")
			{
				colnames(richness) = c("n.taxa", "taxon");
			} else {
				stop("'richness' data appear incorrectly formated: see medusa()")
			}
		}
	}
	
    # checking for typo; if same size, nothing should be dropped
	check <- FALSE;
	if (length(phy$tip.label) == length(richness[,1])) check <- TRUE;
    
    # Prune tree down to lineages with assigned richnesses
	temp <- richness[, "n.taxa"];
	names(temp) <- richness[, "taxon"];
	pruned <- treedata(phy, temp, ...)  # geiger function calling ape (namecheck)
	if (check) {
		if (length(phy$tip.label) != length(pruned$phy$tip.label)) {
			stop("'richness' data and tip labels in 'phy' do not appear to match fully")
		}
	}
	phy <- pruned$phy;
    rr <- pruned$data
    richness <- structure(data.frame(taxon=rownames(rr), n.taxa=rr[,1], row.names=NULL))
    # Check the tree
	#	plotNN(phy)					# Node numbers (ape-style) plotted
    	
	return(list(phy=phy, richness=richness));
}


## Original default was to fit 20 models (or less if the tree was small).
## Changing to a stop-criterion (stop="model.limit") e.g. when k = n-1 (i.e. when denominator of aicc correction is undefined).
## k <- (3*i-1) # when both birth and death are estimated, where i is the number of piecewise models
  ## This occurs when i = n/3
  ## If Yule, max i = n/2
## n <- (2*num.taxa - 1) == (2*length(richness[,1]) - 1) # i.e. total number of nodes in tree (internal + pendant)
## Alternatively use aicc threshold itself as a stopping criterion (stop="threshold").
# AICc = AIC + 2*k*(k+1)/(n-k-1);
.get.max.model.limit <- function (richness, stop, npartitions, model, verbose)
{
	samp.size <- (2*nrow(richness) - 1)
	if (model == "bd" || model == "mixed")
	{
		max.model.limit <- as.integer(samp.size/3) - ((!(samp.size %% 3)) * 1);
	} else {
		max.model.limit <- as.integer(samp.size/2) - ((!(samp.size %% 2)) * 1);
	}
	
	if (stop == "partitions")
	{
		if (npartitions > max.model.limit) {
            model.limit <- max.model.limit;
            warning("Supplied 'partitions' is in excess of the maximal number that can be considered")
        } else {
            model.limit=npartitions
        }
	} else {
		model.limit <- max.model.limit;
	}
	
	if (verbose)
	{
		cat("\nLimiting consideration to a maximum of ", model.limit, " piecewise", sep="");
		if (model == "bd")
		{
			cat(" birth-death models")
		} else if (model == "mixed")
		{
			cat(" mixed models")
		} else {
			cat(" pure-birth (Yule) models");
		}
		if (stop == "threshold") {cat(" (or until threshold is not satisfied)");}
		cat(".\n\n")
	}
	
	return(model.limit);
}


## Fitted curve from random b-d simulations
## Value corresponds to 95th percentile of AICc(split) - AICc(no-split) for no-split simulations
## x-shifted power function
#get.threshold
.threshold.medusa <- function (N)
{
	a = -3.5941052380332650E+01;
	b =  6.7372587299747000E+00;
	c = -1.0061508340754866E-01;
	Offset =  2.7516678664333408E+01;
	y <- a * (N-b)^c + Offset;
	if (y < 0) y <- 0;
	return(y);
}


## The make.cache.medusa function is like the first half of the original splitEdgeMatrix().
## It works through and reorders the edges, then works out start and end times of these
## based on the phylogeny's branching times.
##
## In addition, every node's descendants are also calculated.  The element 'desc' is a list.
## $desc[i] contains the indices within $edge, $t.start, etc., of all descendants of node 'i'
## (in ape node numbering format).
## f: either mclapply or lapply -- see .get.parallel()
.make.cache.medusa <- function (phy, richness, fx)
{
	n.tips <- length(phy$tip.label);
	n.int <- nrow(phy$edge) - n.tips;
	
## Ape numbers the tips first
	i.int <- seq_len(n.int);
	interior <- phy$edge[,2] %in% phy$edge[,1];
	bt <- branching.times(phy);
	
# Consider only internal edges first
	edges.int <- phy$edge[interior,];
	colnames(edges.int) <- c("anc", "dec");
	
	t.0 <- bt[match(edges.int[,1], (n.tips+1):max(edges.int))];
	t.1 <- c(t.0[i.int] - phy$edge.length[interior]);
	
	z.internal <- cbind(edges.int, t.0, t.1, t.len=t.0 - t.1,
		n.0=rep(1, n.int), n.t=rep(NA, n.int));
	
# Now, pendant edges; 
	edges.pendant <- phy$edge[match(seq_len(n.tips), phy$edge[,2]),];
	colnames(edges.pendant) <- c("anc", "dec");
	
	t.0 <- bt[match(edges.pendant[,1], (n.tips+1):max(edges.pendant))];
	t.1 <- rep(0, n.tips);
# cannot assume richness ordering necessarily matches that of tip labels
	ext.richness <- richness$n.taxa[match(phy$tip.label, richness$taxon)];
	
	z.pendant <- cbind(edges.pendant, t.0, t.1, t.len=t.0 - t.1,
		n.0=rep(1, n.tips), n.t=ext.richness);
	
	z <- rbind(z.internal, z.pendant);
	z <- cbind(z,partition=rep(1, length(z[,1]))); # Stores piecewise model structure
	rownames(z) <- NULL;
	
# Used for identifying descendant nodes below i.e. tracking breakpoints
	all.edges <- as.matrix(z[,c("anc","dec")]);
	desc.stem <- list();
	desc.node <- list();

    desc.stem <- fx(seq_len(max(all.edges)), .descendants.cutAtStem.idx, all.edges=all.edges);
    desc.node <- fx(seq_len(max(all.edges)), .descendants.cutAtNode.idx, all.edges=all.edges);

	for (i in 1:length(desc.node))
	{
		if (length(desc.node[[i]]) == 0) # tips
		{
			desc.node[[i]] = .descendants.cutAtStem.idx(node.list=i, all.edges=all.edges)
		}
	}
	
    tips <- .cache.descendants(phy)$tips
	res <- list(z=z, desc.stem=desc.stem, desc.node=desc.node, tips=tips);
	return(res);
}

## This generates the indices of all descendants of a node, using ape's edge matrix.
## Deals with row numbers of the edge matrix rather than node numbers of the tree.
.descendants.cutAtStem <- function (node, all.edges)
{
	ans <- numeric();
	ans <- node;
	repeat {
		node <- all.edges[all.edges[,1] %in% node,2];
		if (length(node) > 0)
		{
			ans <- c(ans, node);
		} else {break;}
	}
	return(unlist(ans));
}

## The function 'descendants' returns the indices of all descendants within the edge matrix.
.descendants.cutAtStem.idx <- function (node.list, all.edges)
{
	which(all.edges[,1] == node.list | all.edges[,2] %in% .descendants.cutAtStem(node.list, all.edges));
}

## This generates the indices of all descendants of a node, using ape's edge matrix.
## Deals with row numbers of the edge matrix rather than node numbers of the tree.
.descendants.cutAtNode <- function (node, all.edges)
{
	ans <- numeric();
	repeat {
		node <- all.edges[all.edges[,1] %in% node,2];
		if (length(node) > 0)
		{
			ans <- c(ans, node);
		} else {break;}
	}
	return(unlist(ans));
}

## The function 'descendants' returns the indices of all descendants within the edge matrix.
.descendants.cutAtNode.idx <- function (node.list, all.edges)
{
	which(all.edges[,1] == node.list | all.edges[,2] %in% .descendants.cutAtNode(node.list, all.edges));
}


## Needed for determining whether nodes are virgin nodes
#get.num.tips
.ntips <- function (node, phy)
{
	n <- length(tips(phy,node));
	return(n);
}


## Only used for base model

#medusa.ml.initial
.initial.medusa <- function (z, initialR, initialE, model)
{
	rootnode <- min(z[,"anc"]);
	obj <- .fit.partition.medusa(partition=1, z=z, sp=c(initialR, initialE), model=model);
	
	model.fit <- .calculate.modelfit.medusa(fit=obj, z=z);
	
	return(list(par=matrix(obj$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik.part=obj$lnLik, 
	   lnLik=obj$lnLik, split.at=rootnode, aic=round(model.fit$aic, digits=7), aicc=round(model.fit$aicc, digits=7), num.par=model.fit$k, cut.at="node"));
}


## Pre-fit values for pendant edges; DON'T recalculate later; should account for ~25% of all calculations
## Also cache values for virgin nodes; useful until subsetted.
## shiftCut can only be "stem" or "node" (not "both"), as both are evaluated separately
#medusa.ml.prefit
.prefit.medusa <- function (node, z, desc, initialR, initialE, model, shiftCut, criterion)
{
	fitted.bd <- NULL;
	fitted.yule <- NULL;
#	fit <- NULL;
	
## if model is mixed, grab the optimally fitted one, drop the other; it is cutAtStem that matters, as it is the sum of 2 break likelihoods
  ## optimal alone may be different than optimal in tandem
	if (shiftCut == "stem")
	{
		if (model == "bd" || model == "mixed")
		{
			obj <- .split.z.at.node.medusa(node=node, z=z, desc=desc, shiftCut="stem");
			z.bd.stem <- obj$z;
			fitted.bd <- .fit.partition.medusa(partition=2, z=z.bd.stem, sp=c(initialR, initialE), model="bd");
			fitted.bd$model <- "bd";
			fitted.bd$cut.at <- "stem";
		}
		if (model == "yule" || model == "mixed")
		{
			obj <- .split.z.at.node.medusa(node=node, z=z, desc=desc, shiftCut="stem");
			z.yule.stem <- obj$z;
			fitted.yule <- .fit.partition.medusa(partition=2, z=z.yule.stem, sp=c(initialR, initialE), model="yule");
			fitted.yule$model <- "yule";
			fitted.yule$cut.at <- "stem";
		}
	} else if (shiftCut == "node")
	{
		if (model == "bd" || model == "mixed")
		{
			obj <- .split.z.at.node.medusa(node=node, z=z, desc=desc, shiftCut="node");
			z.bd.node <- obj$z;
			fitted.bd <- .fit.partition.medusa(partition=2, z=z.bd.node, sp=c(initialR, initialE), model="bd");
			fitted.bd$model <- "bd";
			fitted.bd$cut.at <- "node";
		}
		if (model == "yule" || model == "mixed")
		{
			obj <- .split.z.at.node.medusa(node=node, z=z, desc=desc, shiftCut="node");
			z.yule.node <- obj$z;
			fitted.yule <- .fit.partition.medusa(partition=2, z=z.yule.node, sp=c(initialR, initialE), model="yule");
			fitted.yule$model <- "yule";
			fitted.yule$cut.at <- "node";
		}
	}
## Check which flavour of model fits best
	if (is.null(fitted.bd) & !is.null(fitted.yule))
	{
		return(fitted.yule);
	} else if (is.null(fitted.yule) & !is.null(fitted.bd)) {
		return(fitted.bd);
	} else {
## Dealing with a 'mixed' model here; need to consider number of parameters.
		bd.model.fit <- .calculate.modelfit.medusa(fit=fitted.bd, z=z);
		yule.model.fit <- .calculate.modelfit.medusa(fit=fitted.yule, z=z);
		element <- NULL;
		if (criterion == "aic") {element <- 1;} else {element <- 2;}
		if (is.nan(yule.model.fit[[element]]))
		{
			return(fitted.bd);
		} else if (bd.model.fit[[element]] < yule.model.fit[[element]])
		{
			return(fitted.bd);
		} else {
			return(fitted.yule);
		}
	}
}


## sp = initializing values for r & epsilon
## Default values should never be used (except for first model), as the values from the previous model are passed in
#medusa.ml.fit.partition
.fit.partition.medusa <- function (partition, z, sp=c(0.1, 0.05), model)
{
# Construct likelihood function:
	lik <- .lik.partition.medusa(partition=(z[z[,"partition"] == partition,,drop=FALSE]), model=model);
	foo <- function (x) {-lik(pars=exp(x));} # work with parameters in log-space to preserve precision
	
	if (model == "bd")
	{
		fit <- optim(fn=foo, par=log(sp), method="N"); # last argument connotes maximization
		return(list(par=exp(fit$par), lnLik=-fit$value));
	} else {
		fit <- optimize(f=foo, interval=c(-25, 1));
		par <- c(exp(fit$minimum), NA);
		return(list(par=par, lnLik=-fit$objective));
	}
}


## Split the edge matrix 'z' by adding a partition rooted at node 'node'.
##   Note: in original MEDUSA parlance, this is cutAtStem=T.
## The list 'desc' is a list of descendants (see make.cache.medusa, above).
## Returns a list with elements:
##   z: new medusa matrix, with the new partition added
##   affected: indices of the partitions affected by the split (n == 2).
## This is where 'shiftCut' matters
#medusa.split
.split.z.at.node.medusa <- function (node, z, desc, shiftCut)
{
	descendants <- NULL;
	if (shiftCut == "stem")
	{
		descendants <- desc$desc.stem;
	} else {
		descendants <- desc$desc.node;
	}
	part <- z[,"partition"];
	base <- min(part[z[,1] == node | z[,2] == node]);
	tag <- max(part) + 1;
	i <- descendants[[node]];
	idx <- i[part[i] == base];
	z[idx,"partition"] <- tag;
	z[which(z["dec"] == node),"partition"] <- tag; # Possible to have several edges to consider
	return(list(z=z, affected=c(unique(part[idx]), tag)));
}


## 'fit' contains parameter values from previous model, used to initialize subsequent model.
## Pass in pre-fitted values for pendant edges and virgin nodes (in 'prefit'); DON'T recalculate.
## Need to consider the possibility of birth-death, yule, or mixed models.
## Need to consider where shft is placed (shiftCut). Placement affects not only new clade, but
## also the size of the split clade. Only relevant if shiftCut = "both".
## fit1 model is already logged; only need to record fit2 model, and only non-prefitted nodes
#medusa.ml.update
.update.fit.medusa <- function (node, z, desc, fit, prefit, num.tips, root.node, model, criterion, shiftCut)
{
## various combinations possible
	fit1.stem <- NULL;
	fit1.node <- NULL;
	fit2.stem <- NULL;
	fit2.node <- NULL;
	cut.at <- NULL;
	
	sp <- NULL;
	aff <- NULL;
	op <- fit$par;
	cut.at <- NULL;
	
	fit1 <- NULL;
	fit2 <- NULL;
	
	if (shiftCut == "stem" || shiftCut == "both")
	{
## First, diminshed clade
		obj.stem <- .split.z.at.node.medusa(node=node, z=z, desc=desc, shiftCut="stem");
		z.stem <- obj.stem$z;
		aff <- obj.stem$affected;
		sp <- op[aff[1],]; # Use previously fit parameter values from clade that is currently being split
		
		if (model == "mixed") ## In mixed models, want to conserve flavour of previously fit model (right?)
		{
			if (sum(!is.na(sp)) < 2) # yule
			{
				fit1.stem <- .fit.partition.medusa(partition=aff[1], z=z.stem, sp=sp, model="yule");
			} else {
				fit1.stem <- .fit.partition.medusa(partition=aff[1], z=z.stem, sp=sp, model="bd");
			}
		} else {
			fit1.stem <- .fit.partition.medusa(partition=aff[1], z=z.stem, sp=sp, model=model);
			fit1.stem$model <- model;
		}
## Second, new clade
		if (node < root.node) # tip, already calculated
		{
			fit2.stem <- prefit$tips[[node]];
		} else if (length(unique(z.stem[(z.stem[,"partition"] == aff[2] & z.stem[,"dec"] < root.node),"dec"])) == num.tips[[node]]) {
			fit2.stem <- prefit$virgin.nodes$stem[[node - root.node]];
		} else { # novel shift
			fit2.stem.bd <- NULL;
			fit2.stem.yule <- NULL;
			
			if (model == "yule" || model == "mixed")
			{
				fit2.stem.yule <- .fit.partition.medusa(aff[2], z.stem, sp=sp, model="yule");
				fit2.stem.yule$model <- "yule";
			}
			if (model == "bd" || model == "mixed")
			{
				if (is.na(sp[2])) {sp[2] <- 0.5;}
				fit2.stem.bd <- .fit.partition.medusa(aff[2], z.stem, sp=sp, model="bd");
				fit2.stem.bd$model <- "bd";
			}
## Figure out which model fits best
			if (is.null(fit2.stem.bd))
			{
				fit2.stem <- fit2.stem.yule;
			} else if (is.null(fit2.stem.yule)) {
				fit2.stem <- fit2.stem.bd;
			} else {
## Considering both places for a shift
				fit2.stem.bd.val <- .calculate.modelfit.medusa(fit=fit2.stem.bd, z=z);
				fit2.stem.yule.val <- .calculate.modelfit.medusa(fit=fit2.stem.yule, z=z);
				
				if (criterion == "aic") {element <- 1;} else {element <- 2;}
				
				if (fit2.stem.bd.val[[element]] < fit2.stem.yule.val[[element]])
				{
					fit2.stem <- fit2.stem.bd;
				} else {
					fit2.stem <- fit2.stem.yule;
				}
			}
		}
	}
	if (shiftCut == "node" || shiftCut == "both")
	{
## First, diminshed clade
		obj.node <- .split.z.at.node.medusa(node=node, z=z, desc=desc, shiftCut="node");
		z.node <- obj.node$z;
		aff <- obj.node$affected;
		sp <- op[aff[1],]; # Use previously fit parameter values from clade that is currently being split
		
		if (model == "mixed") ## In mixed models, want to conserve flavour of previously fit model (right?)
		{
			if (sum(!is.na(sp)) < 2) # yule
			{
				fit1.node <- .fit.partition.medusa(partition=aff[1], z=z.node, sp=sp, model="yule");
				fit1.node$model <- "yule";
			} else {
				fit1.node <- .fit.partition.medusa(partition=aff[1], z=z.node, sp=sp, model="bd");
				fit1.node$model <- "bd";
			}
		} else {
			fit1.node <- .fit.partition.medusa(partition=aff[1], z=z.node, sp=sp, model=model);
			fit1.node$model <- model;
		}
## Second, new clade
		if (node < root.node) # tip, already calculated
		{
			fit2.node <- prefit$tips[[node]];
		}
		else if (length(unique(z.node[(z.node[,"partition"] == aff[2] & z.node[,"dec"] < root.node),"dec"])) == num.tips[[node]]) # virgin node
		{
			fit2.node <- prefit$virgin.nodes$node[[node - root.node]];
		}
		else { # novel shift
			fit2.node.bd <- NULL;
			fit2.node.yule <- NULL;
			
			if (model == "yule" || model == "mixed")
			{
				fit2.node.yule <- .fit.partition.medusa(aff[2], z.node, sp, model="yule");
				fit2.node.yule$model <- "yule";
			}
			if (model == "bd" || model == "mixed")
			{
				if (is.na(sp[2])) {sp[2] <- 0.5;}
				fit2.node.bd <- .fit.partition.medusa(aff[2], z.node, sp, model="bd");
				fit2.node.bd$model <- "bd";
			}
## Figure out which model fits best
			if (is.null(fit2.node.bd))
			{
				fit2.node <- fit2.node.yule;
			} else if (is.null(fit2.node.yule)) {
				fit2.node <- fit2.node.bd;
			} else {
## Considering both places for a shift
				fit2.node.bd.val <- .calculate.modelfit.medusa(fit=fit2.node.bd, z=z);
				fit2.node.yule.val <- .calculate.modelfit.medusa(fit=fit2.node.yule, z=z);
				
				if (criterion == "aic") {element <- 1;} else {element <- 2;}
				
				if (fit2.node.bd.val[[element]] < fit2.node.yule.val[[element]])
				{
					fit2.node <- fit2.node.bd;
				} else {
					fit2.node <- fit2.node.yule;
				}
			}	
		}
	}
	
## Now, figure out which shift position is optimal	
	if (is.null(fit2.node))
	{
		fit1 <- fit1.stem;
		fit2 <- fit2.stem;
		cut.at <- "stem";
	} else if (is.null(fit1.stem)) {
		fit1 <- fit1.node;
		fit2 <- fit2.node;
		cut.at <- "node";
	} else {
## Considering both places for a shift
		stem.lik <- (fit1.stem$lnLik + fit2.stem$lnLik);
		stem.par <- rbind(fit1.stem$par, fit2.stem$par)
		stem.val <- list(lnLik=stem.lik, par=stem.par);
		stem.fit <- .calculate.modelfit.medusa(fit=stem.val, z=z);
		
		node.lik <- (fit1.node$lnLik + fit2.node$lnLik);
		node.par <- rbind(fit1.node$par, fit2.node$par)
		node.val <- list(lnLik=node.lik, par=node.par);
		node.fit <- .calculate.modelfit.medusa(fit=node.val, z=z);
		
#		element <- NULL;
#		if (criterion == "aic") {element <- 1;} else {element <- 2;}
		
		if (stem.fit[[criterion]] < node.fit[[criterion]])
		{
			fit1 <- fit1.stem;
			fit2 <- fit2.stem;
			cut.at <- "stem";
#			cat("Stem wins!\n");
		} else {
			fit1 <- fit1.node;
			fit2 <- fit2.node;
			cut.at <- "node";
#			cat("Node wins!\n");
		}
	}
	op[aff[1],] <- fit1$par; # Replace parameters with new values for diminished clade
	
	fit$par <- rbind(op, fit2$par);
	fit$lnLik.part[aff] <- c(fit1$lnLik, fit2$lnLik); # Replace parameters with new values for diminished clade
	fit$split.at <- c(fit$split.at, node);
	fit$lnLik <- sum(fit$lnLik.part);
	
	model.fit <- .calculate.modelfit.medusa(fit=fit, z=z);
	
	fit$aic <- model.fit$aic;
	fit$aicc <- model.fit$aicc;
	fit$num.par <- model.fit$k;
	fit$cut.at <- cut.at;
	fit$model <- fit2$model;
	
	return(fit);
}


## make.lik.medusa.part: generate a likelihood function for a single partition.
#make.lik.medusa.part
.lik.partition.medusa <- function (partition, model)
{
# Handle internal and pendant edges separately
	is.int <- is.na(partition[,"n.t"]);
	is.pend <- !is.int;
	
	n.int <- sum(is.int);
	n.pend <- sum(is.pend);
	
	if (n.int + n.pend != length(partition[,1])) stop("You messed up, yo.");
	
## Internal and pendant calculations differ; split'em up
	int  <- partition[is.int,,drop=FALSE];
	pend <- partition[is.pend,,drop=FALSE];
	
	sum.int.t.len <- sum(int[,"t.len"]);  # Simply sum all internal edges
	int.t.0 <- int[,"t.0"];
	
# 'n.0' = Foote's 'a', initial diversity; 'n.t' = Foote's 'n', final diversity
	pend.n.0 <- pend[,"n.0"]; # Foote's 'a': initial diversity
	pend.n.t <- pend[,"n.t"]; # Foote's 'n': final diversity
	pend.t.len <- pend[,"t.len"];
	
# User may pass in epsilon; don't change it, just estimate r
	f <- function(pars)
	{
		if (model == "bd")
		{
			r <- pars[1];
			epsilon <- pars[2];
			
			if (r < 0 || epsilon <= 0 || epsilon >= 1) {return(-Inf);}
		} else {
			
			r <- pars[1];
			epsilon <- 0;
			
			if (r < 0) {return(-Inf);}
		}
		
		l.int <- numeric();
		l.pend <- numeric();
		
		if (n.int == 0) {l.int <- 0;} else {
## Likelihood of internal edges from Rabosky et al. (2007) equation (2.3):
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
		}
		
		if (n.pend == 0) {l.pend <- 0} else {
## Calculations are from the following:
## Rabosky et al. 2007. Proc. Roy. Soc. 274: 2915-2923.
## Foote et al. 1999. Science. 283: 1310-1314
## Raup. 1985. Paleobiology 11: 42-52 [Foote et al. correct the equation [A18] where a > 1]
## Bailey. 1964. The Elements Of Stochastic Processes, With Applications To The Natural Sciences
## Kendall. 1948. Ann. Math. Stat. 19: 1–15.
##
## A = probability of extinction of one lineage over time 't'
## B = A * (lambda/mu)
##
## When there is a single lineage at time 0 (a = 1), the calculation is
##   log(1 - A) + log(1 - B) + (n - 1)*log(B)
## but this is conditioned on survival by dividing by (1-A)
## (subtracting log(1 - A) on a log scale) which cancels to give:
##   log(1 - B) + (n - 1)*log(B)
##      - for n.t == 1, reduces further to log(1-B)
##
## A = mu*(exp((lambda - mu)*t) - 1)) / (lambda*exp((lambda - mu)*t) - mu)
##  let r = (lambda - mu); ert = exp((lambda - mu)*t)
## A = mu*(ert - 1)/(lambda*ert - mu)
##
## B = A * (lambda/mu)
##   = [mu*(ert - 1)/(lambda*ert - mu)] * (lambda/mu)
##   = (lambda*(ert - 1))/(lambda*ert - mu)
##   = (lambda*(ert - 1))/(lambda(ert - mu/lambda))
##   = (ert - 1) / (ert - epsilon)

## All pendant nodes begin with richness '1'; calculations simple.
#			i.pend.n.t.1 <- which(pend.n.t == 1)   # calculations even simpler: log(1-B)
#			i.pend.n.t.n1 <- which(pend.n.t != 1)
			
			ert <- exp(r * pend.t.len);
			B <- (ert - 1) / (ert - epsilon); # Equivalently: B <- (bert - b) / (bert - d)
			
			l.pend <- sum(log(1 - B) + (pend.n.t - 1)*log(B));
		}
		return(l.int + l.pend);
	}
}


## 'fit' contains '$par' and '$lnlik'
#calculate.model.fit
.calculate.modelfit.medusa <- function (fit, z)
{
## Sample size taken (for now) as the total num.nodes in the tree (internal + pendant)
  # num.nodes = (2*length(phy$tip.label) - 1) == (2*length(richness[,1]) - 1) == length(z[,1]) + 1
#	n <- (length(z[,1]) + 1) + sum(!is.na(z[,"n.f"]));
	
# Since each edge defines a node (i.e. an 'observation'), need only add root node as final obervation
	n <- (length(z[,1]) + 1);
	
 # Includes both formal parameters AND number of breaks. Note: first model does not involve a break.
## Models where all parameters are estimated (i.e. BD model):
  # 2 parameters for base model (no breakpoint) + 3 parameters (r, eps, breakpoint) for each subsequent model
  
  
# Determine number of piecewise models currently involved
	if (length(fit$par) < 3) # i.e. base model
	{
		num.models <- 1;
	} else {
		num.models <- length(fit$par[,1]);
	}
	
# Updated for more general models: check how many parameter values != NA
#	k <- 2 + (3 * (num.models - 1))
	k <- sum(!is.na(fit$par)) + (num.models - 1); # number of estimated parameters + number of breaks
	
    ll = list(k=k, lnL=fit$lnLik)
    return(.aic(ll, n))
    #	lnLik <- fit$lnLik;
	
    #	aic <- (-2 * lnLik) + (2*k);
    #	aicc <- aic + 2*k*(k+1)/(n-k-1);
	
    #	model.fit <- c(aic, aicc, k);
    #	return(model.fit);
}


## Prints out a table of likelihoods, parameters, and aic scores
# calculate.model.fit.summary
.summary.modelfit.medusa <- function (models, phy, threshold=0)
{
	tmp <- matrix(nrow=(length(models)), ncol=6);
	colnames(tmp) <- c("partitions", "split", "lnL", "k", "aic", "aicc");
    # c("N.Models", "Break.Node", "Ln.Lik", "N.Param", "aic", "aicc")
	
	w.aic <- numeric(length(models));
	w.aicc <- numeric(length(models));
	cut.at <- character(length(models));
	
	for (i in 1:length(tmp[,1]))
	{
		tmp[i,] <- c(i, as.integer(models[[i]]$split.at[i]), models[[i]]$lnLik, models[[i]]$num.par, models[[i]]$aic, models[[i]]$aicc);
		cut.at[i] <- models[[i]]$cut.at;
	}
	
	all.res <- data.frame(tmp);
	all.res[1,2] <- NA # root node for base model
    all.res$cut=as.character(cut.at)
    all.res=all.res[,c("partitions", "split", "cut", "lnL", "k", "aic", "aicc")]

	if (threshold == 0)
	{
		w.aic <- round(.daic(all.res$aic), digits=5);
		w.aicc <- round(.daic(all.res$aicc), digits=5);
        all.res$w.aic=w.aic$w
        all.res$w.aicc=w.aicc$w
 	} 
	
    #	if (plot)
    #	{
        #		dev.new();
    #		.plot.modelfit.medusa(all.res, ...);
        #		if (!is.null(fig.title)) {title(main=fig.title, cex.main=0.75);}
    #	}
	return(all.res);
}


## Self explanatory
## These are meaningless when using a threshold criterion
#calculate.model.weights
.daic <- function (aic)
{
	best <- min(aic);
	delta <- aic-best;
	sumDelta <- sum(exp(-0.5 * delta));
	w <- (exp(-0.5 * delta)/sumDelta);
	
	results <- data.frame(fit=aic,delta=delta,w=w);
    rownames(results)=names(aic)
	
	return(results);
}

summary.medusaRAW <- function(object, criterion=c("aicc", "aic"), model=NULL, ...)
{
    #    function(results, modelNum=NULL, cutoff="threshold", criterion="aicc", plotTree=TRUE, time=TRUE, node.labels=TRUE, cex=0.5, plotSurface=FALSE, n.points=100, ...)
    # Desirables:
    #  1. table listing parameter values of selected model
    #  2. list parameters of base model
    #  3. tree printed with colour-coded edges, node labels to indicate split position(s)
    #  4. plot likelihood surface
	
    # Extract constituent components from results
    #    ct=list(model=NULL, criterion=c(NA, "aic", "aicc"), threshold=0)
    #    ct[names(control)]=control
    
    criterion=match.arg(criterion, c("aicc", "aic"))
    
    results <- object
	fit <- results$models;
	phy <- results$cache$phy;
	desc <- results$cache$desc;
	modelSummary <- results$summary;
    cutAt <- as.character(modelSummary$cut);
	
    #	fit; phy; z; desc; modelSummary; threshold;
	
    # First, determine which model is desired
    
    model.chooser=function(criterion=c("aic", "aicc"), threshold="auto"){
        if(threshold=="auto"){
            threshold=.threshold.medusa(Ntip(phy))
        } else if(!is.numeric(threshold)) {
            stop("'threshold' must be numeric")
        }
        threshold=abs(threshold)
        
        model.id <- 1;
		while (1)
		{
			if ((model.id + 1) > length(fit)) break;
			if ((fit[[model.id]][[criterion]] - fit[[model.id+1]][[criterion]]) < threshold) break;
			model.id <- model.id + 1;
		}
        structure(threshold, names=model.id)
    }
    
    
	model.id <- 0;
	if (!is.null(model))
	{
        xx=list(...)
        if("threshold"%in%names(xx)) warning("'threshold' has been ignored")
		model.id <- model;
        threshold = 0
	} else {   # Find best model using some criterion (threshold or user-defined)
        #		if (cutoff != "threshold") {threshold <- cutoff}
        #		else {cat("\nSelecting model based on corrected threshold (improvement in information theoretic score of ", ct$threshold, " units).\n", sep="");}
        tmp=model.chooser(criterion=criterion, ...)
        model.id=as.integer(names(tmp))
        threshold=as.numeric(tmp)
        
	}
    
    names(threshold)=criterion
    chosen=ifelse(1:length(fit)%in%model.id, "*", "")
    modelSummary$chosen=chosen
    z=object$FUN(model.id)
    attr(z, "phylo")=phy
    attr(z, "summary")=modelSummary
    attr(z, "threshold")=threshold
    class(z)=c("medusa", class(z))
    return(z)
}


plot.medusa = function(x, partitions=list(cex=2, bg="gray", alpha=0.75, col="black", lwd=1), ...){
    z=x
    st=list(cex=2, bg="gray", alpha=0.75, col="black", lwd=1)
    st[names(partitions)]=partitions
    shift=st
    phy=attr(z, "phylo")
    #    par=match.arg(par, c("r", "epsilon"))
    #    bb=matrix(z[,par], nrow=1)
    #    bb=rbind(bb, bb)
    #    colnames(bb)=rownames(z)
    #    shifts=max(z[,"partition"])
    #    yy=.branchcol.plot(phy, bb, plot=FALSE, colors=list(branches=min(c(256, 4*shifts+1)), legend=7, missing=1))$col
    plot(phy, ...)
    shifts=z[(idx<-!is.na(z[,"shift"])),"dec"]
    if(length(shifts)){
        ss=z[idx,"shift"]
        ww=character(nrow(phy$edge))
        ww[match(shifts, phy$edge[,2])]=ss
        xx=numeric(nrow(phy$edge))
        xx[phy$edge[,2]%in%shifts]=1
        edgelabels.auteur(NULL, frame="circle", cex=ifelse(xx==1, shift$cex, 1e-10), pch=ifelse(xx==1, 21, NA), bg=.transparency(shift$bg, shift$alpha), col=shift$col, lwd=shift$lwd)
        edgelabels.auteur(ww, frame="none", cex=ifelse(xx==1, shift$cex/3, 1e-10), pch=NA, col=shift$col)
    }
}

plot.medusaRAW = function(x, col=c(aic="black", aicc="blue"), ...){
    all.res=x$summary
    ylim <- c(min(all.res[,"aic"],all.res[,"aicc"]), max(all.res[,"aic"],all.res[,"aicc"]));
	plot(all.res[,"partitions"], all.res[,"aicc"], xlab="piecewise models", ylab="model fit", ylim=ylim, type="l", col=col[["aicc"]], xaxt="n", ...);
    axis(1, at=all.res[,"partitions"], labels=all.res[,"partitions"], tcl=NA)
	points(all.res[,"partitions"],all.res[,"aicc"], col=col[["aicc"]], pch=21, bg="white");
	points(all.res[,"partitions"],all.res[,"aic"], col=col[["aic"]], type="l");
	points(all.res[,"partitions"],all.res[,"aic"], col=col[["aic"]], pch=21, bg="white");
	
	legend("topleft", c("aicc","aic"), pch=21, pt.bg="white", lty=1, col=c(col[["aicc"]], col[["aic"]]), inset = .05, cex=0.75, bty="n");
}



## Create a plot of model-fit vs. model-size
#plotModelFit
#.plot.modelfit.medusa <- function (all.res, ...)
#{
#	ylim <- c(min(all.res[,"aic"],all.res[,"aicc"]), max(all.res[,"aic"],all.res[,"aicc"]));
#	plot(all.res[,"partitions"],all.res[,"aicc"], xlab="number of piecewise models", ylab="model fit", ylim=ylim, type="l", col="blue", xaxt="n", ...);
#   axis(1, at=all.res[,"partitions"], labels=all.res[,"partitions"], tcl=NA)
#	points(all.res[,"partitions"],all.res[,"aicc"], col="blue", pch=21, bg="white");
#	points(all.res[,"partitions"],all.res[,"aic"], col="black", type="l");
#	points(all.res[,"partitions"],all.res[,"aic"], col="black", pch=21, bg="white");
	
#	legend("topleft", c("aicc","aic"), pch=21, pt.bg="white", lty=1, col=c("blue", "black"), inset = .05, cex=0.75, bty="n"); # 'bottomright' also works
#}

## Takes in a summary from summarizeTurboMEDUSA
## treeParameters <- list(mm=mm, break.pts=break.pts, phy=phy, z=z)
#plotPrettyTree <- function (treeParameters, time=TRUE, node.labels=FALSE, margin=FALSE, cex=0.5, label.offset=0, font=3, color.tip.label=FALSE, ...)
#{
#	mm <- treeParameters$mm;
#	break.pts <- treeParameters$break.pts;
#	phy <- treeParameters$phy;
#	z <- treeParameters$z;
#	colour <- NULL;
	
#	dev.new();
#
#	if (color.tip.label)
#	{
#		for (i in 1:length(phy$tip.label))
#		{
#			colour[i] <- as.integer(z[which(z[,"dec"] == i),"partition"]);
#		}
#	}
#	if (time) {margin=TRUE;}
#	plot.phylo(phy, edge.color=z[mm,"partition"], no.margin=!margin, cex=cex, label.offset=label.offset, tip.color=colour, ...);
#	if (time)
#	{
#		axisPhylo(cex.axis=0.75);
#		mtext("Divergence Time (MYA)", at=(max(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx)*0.5), side = 1, line = 2, cex=0.75);
#	}
#	if (node.labels)
#	{
#		for (i in  1:length(break.pts))
#		{
#			nodelabels(i, node= break.pts[i], frame = "c", font = 1, cex=0.5);
#		}
#	}
#}


## Get b and d values from r (b-d) and epsilson (d/b)
## Used in previous version of program; now in terms of r and epsilon
## Possibly of use to users wishing to translate results
#get.b.d <- function (r, epsilon)
#{
#	b <- r/(1-epsilon);
#	d <- b-r;   # Alternatively: d <- eps*r/(1-eps)
#	return(list(b=b, d=d));
#}

## Print out tree with ape-style node-numbering
## Possibly of interest for users to identify numbers of node(s) off interest
 ## If this is the case, make sure to pass in pruned tree
#plotNN <- function (phy, time=TRUE, margin=TRUE, label.offset=0.5, cex=0.5, ...)
#{
#	phy$node.label <- (length(phy$tip.label) + 1):max(phy$edge);
#	plot.phylo(phy, show.node.label=TRUE, no.margin=!margin, label.offset=label.offset, cex=cex, ...);
#	if (time && !margin) {cat("Cannot plot time axis without a margin.\n");}
#	else if (time && margin) {axisPhylo(cex.axis=0.75)};
#}