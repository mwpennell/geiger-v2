medusaVersion = 1.05;

medusa <- function (phy, richness = NULL, criterion = c("aicc", "aic"), partitions = NA, model = c("mixed", "bd", "yule"),
	cut = c("both", "stem", "node"), stepBack = TRUE, init = c(r = 0.05, epsilon = 0.5), ncores = NULL, verbose = FALSE, ...) {

	## CHECK ARGUMENTS
	#verbose <- FALSE; # this should be an argument above
	initialE <- init[["epsilon"]];
	initialR <- init[["r"]];
	sp = c(initialR, initialE);

	fx <- .get.parallel(ncores);

	if (!is.na(partitions)) {
		flag <- "'partitions' should either be NA or an integer, specifying the maximum number of piecewise models to consider";
		if (!is.numeric(partitions)) {
			stop(flag);
		}
		if (partitions <= 0) {
			stop(flag);
		}
		stop <- "partitions";
	} else {
		criterion <- match.arg(criterion, choices = c("aicc", "aic"), several.ok = FALSE);
		stop <- "threshold";
	}
	if (stop == "threshold") {
		npartitions <- 0;
	} else {
		npartitions <- partitions;
		criterion <- "aicc";
	}

	model <- match.arg(model, choices = c("mixed", "bd", "yule"), several.ok = FALSE);
	shiftCut <- match.arg(cut, choices = c("both", "stem", "node"), several.ok = FALSE);
	
	## Before determining model.limit, prune tree as necessary (from 'taxon' information in 'richness')
	if (!any(c("phylo", "multiPhylo") %in% class(phy))) {
		stop("'phy' must either be a phylo or multiPhylo object");
	}
	richness <- .check.richness(phy = phy, richness = richness);
	phyData <- .treedata.medusa(phy = phy, richness = richness, warnings = FALSE); ## modified prune.tree.merge.data for multiple trees (jme)
	
	# set threshold/model.limit once instead of for every tree
	threshold_N <- ifelse(stop == "threshold", .threshold.medusa(phyData$phy), 0);
	model.limit <- .get.max.model.limit(richness = richness, stop = stop, npartitions = npartitions, model = model, verbose = verbose);
	
	# internal function for running a single tree and richness data.frame (jme)
	medusa_runner <- function (phy, richness) {
		
	## move both of the following 2 outside of the function. arg.
		
		## Limit on number of piecewise models fitted; based on tree size, aicc correction factor,
		## and flavour of model fitted (i.e. # parameters estimated; birth-death or pure-birth)
		#model.limit <- .get.max.model.limit(richness = richness, stop = stop, npartitions = npartitions, model = model, verbose = verbose);
		
		## Determine correct AICc threshold from tree size (based on simulations)
		## Should be used for interpreting model-fit
		#N <- Ntip(phy);
		#threshold_N <- ifelse(stop == "threshold", .threshold.medusa(N), 0);
		#cat("Appropriate AICc threshold for tree of ", N, " tips is: ", threshold_N, ".\n\n", sep="");
		
		## Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
		pend.nodes <- seq_len(length(phy$tip.label)); # Calculate pendant splits just once, keep track through various models
		int.nodes <- unique(phy$edge[,1])[-1]; # Omit root node
		root.node <- length(phy$tip.label) + 1;
		all.nodes <- c(pend.nodes, root.node, int.nodes);
		
		## Store pertinent information: branch times, richness, descendants
		#cat("Preparing data for analysis... ");
		obj <- .make.cache.medusa(phy = phy, richness = richness, fx = fx, shiftCut = shiftCut);
		#cat("done.\n");
		desc <- list(stem = obj$desc.stem, node = obj$desc.node);
		
		# update to this:
		#desc <- list(stem = obj$desc.stem, node = obj$desc.node);
		
		# WTF is z.orig doing here? motherfucker...
		z <- z.orig <- obj$z;
		
		## Fit the base model
		## 'fit' holds current results; useful for initializing subsequent models
		
		
		## update this to use medusaFitOptimal
		fit <- list();
		fit <- .fit.base.medusa (z = z, sp = sp, model = model, criterion = criterion);
		
		# if (model == "mixed") {
			# fit.bd <- .base.medusa(z = z, sp = sp, model = "bd");
			# fit.yule <- .base.medusa(z = z, sp = sp, model = "yule");
			# if (fit.bd[[criterion]] < fit.yule[[criterion]]) {
				# fit <- fit.bd;
				# fit$model <- "bd";
			# } else {
				# fit <- fit.yule;
				# fit$model <- "yule";
			# }
		# } else if (model == "bd") {
			# fit <- .base.medusa(z = z, sp = sp, model = "bd");
			# fit$model <- "bd";
		# } else {
			# fit <- .base.medusa(z = z, sp = sp, model = "yule");
			# fit$model <- "yule";
		# }
		#fit$z <- z;

		#models <- list(fit); # get rid of this bullshit...
		#zz <- list(z); # get rid of this bullshit...
		
	# If only one model is desired (i.e. base model), don't bother with all of the precalculations.
		prefit <- NULL;
		if (model.limit != 1) {
			
		## Needed downstream; do not recalculate
		## Gives the number of tips associated with an internal node; determines whether a node is 'virgin' or not
			num.tips <- list();
			num.tips <- fx(all.nodes, function (x) length(obj$tips[[x]]));
		
		## Pre-fit pendant edges so these values need not be re(re(re))calculated; amounts to ~25% of all calculations
		## Will show particular performance gain for edges with many fossil observations
			#cat("Optimizing parameters for pendant edges... ");
			tips <- NULL;
			# Will always be shiftCut="stem"; if mixed model, keep only best fit and throw out other in .prefit.medusa
			#tips <- fx(pend.nodes, .prefit.medusa, z = z, desc = desc, sp = sp, model = model, shiftCut = "stem", criterion = criterion);
			tips <- fx(pend.nodes, .prefit.tip.medusa, z = z, sp = sp, model = model, criterion = criterion);
			#cat("done.\n");
			
		## Pre-fit virgin internal nodes; should deliver performance gain for early models, and especially for large trees
		## Remain useful until a split is accepted within the clade
			## Need to incorporate cutAtStem here
			#cat("Pre-calculating parameters for virgin internal nodes... ");
			virgin.stem <- list();
			virgin.node <- list();
	
			if (shiftCut == "stem" || shiftCut == "both") {
				#virgin.stem <- fx(int.nodes, .prefit.medusa, z = z, desc = desc, sp = sp,
				#	model = model, shiftCut = "stem", criterion = criterion);
				virgin.stem <- fx(int.nodes, .prefit.medusa, z = z, desc = desc, sp = sp,
					model = model, shiftCut = "stem", criterion = criterion);
			}
			if (shiftCut == "node" || shiftCut == "both") {
				#virgin.node <- fx(int.nodes, .prefit.medusa, z = z, desc = desc, sp = sp,
				#	model = model, shiftCut = "node", criterion = criterion);
				virgin.node <- fx(int.nodes, .prefit.medusa, z = z, desc = desc, sp = sp,
					model = model, shiftCut = "node", criterion = criterion);
			}
			
			virgin.nodes <- list(stem = virgin.stem, node = virgin.node);
			#cat("done.\n\n");
			prefit <- list(tips = tips, virgin.nodes = virgin.nodes);
		}
		
		# this is turned off for the moment
		#if (stop == "partitions") {
		stopOnLimit <- function(fit) { # not likely to be used. get rid of it?
			optModel <- fit;
			i <- 1;
			
			#cat("Step 1 (of ", model.limit, "): best likelihood = ", models[[1]]$lnLik, "; AICc = ", models[[1]]$aicc,
			#	"; model = ", models[[1]]$model, "\n", sep="");
			cat("Step 1: lnLik=", round(optModel$lnLik, digits=7), "; ", criterion, "=",
				round(as.numeric(optModel[criterion]), digits=7), "; model=", optModel$model[1], "\n", sep="");
			
			#for (i in seq_len(model.limit - 1)) {
			while (i < model.limit) {
				i <- i + 1;
				node.list <- all.nodes[-fit$split.at];
				res <- fx(node.list, .update.fit.medusa, z = z, desc = desc, fit = fit, prefit = prefit, num.tips = num.tips,
					root.node = root.node, model = model, criterion = criterion, shiftCut = shiftCut);

			# Select model with best score according to the specific criterion employed (default aicc)
				#best <- which.min(unlist(lapply(res, "[[", criterion)));
				fit <- res[[which.min(unlist(lapply(res, "[[", criterion)))]];
				
				#models <- c(models, res[best]);
				#fit <- res[[best]]; # keep track of '$split.at' i.e. nodes already considered
				#z <- zz[[i + 1]] <- .split.z.at.node.medusa(node = node.list[best], z = z, desc = desc, shiftCut = fit$cut.at)$z;
				
				z.best <- .split.z.at.node.medusa(node = tail(fit$split.at, 1), z = z, desc = desc, shiftCut = tail(fit$cut.at, 1))$z;
				step <- rbind(optModel$step, c("add", tail(fit$split.at, 1)));
				
			# Consider parameter removal
				if (stepBack) {
					backFit <- .back.step.medusa(currentModel = fit, z = z.best, step = step, model = model, criterion = criterion);
					fit <- backFit$fit;
					z.best <- backFit$z;
					step <- backFit$step;
					if (!is.null(backFit$remove)) {
						nn <- length(backFit$remove);
						i <- i - nn;
					}
				}
				fit$step <- step;
				
				if (stepBack && !is.null(backFit$remove)) {.print.removed.shifts(remove=backFit$remove);}
				optModel <- fit;
				z <- z.best;
				
				#cat("Step ", i+1, " (of ", model.limit, "): best likelihood = ", round(models[[i+1]]$lnLik, digits=7),
				#	"; AICc = ", models[[i+1]]$aicc, "; break at node ", models[[i+1]]$split.at[i+1], "; model=", models[[i+1]]$model,
				#	"; cut=", models[[i+1]]$cut.at, "\n", sep="");
				cat("Step ", i, ": lnLik=", round(fit$lnLik, digits=7),
					"; AICc=", fit$aicc, "; shift at node ", tail(fit$split.at,1), "; model=",
					tail(fit$model,1), "; cut=", tail(fit$cut.at,1), "; # shifts=", length(fit$split.at) - 1, "\n", sep="");
			}
			return(optModel);
		}
		
		#else if (stop == "threshold") {
		stopOnThreshold <- function (fit) {
			optModel <- fit;
			i <- 1;
			done <- FALSE;
			
			#cat("Step 1: best likelihood = ", models[[1]]$lnLik, "; AICc = ", models[[1]]$aicc, "; model = ", models[[1]]$model, "\n", sep="");
			cat("Step 1: lnLik=", round(optModel$lnLik, digits=7), "; ", criterion, "=",
				round(as.numeric(optModel[criterion]), digits=7), "; model=", optModel$model[1], "\n", sep="");
			
			while (!done & i < model.limit) {
				node.list <- all.nodes[-fit$split.at];
				res <- fx(node.list, .update.fit.medusa, z = z, desc = desc, fit = fit, prefit = prefit, num.tips = num.tips,
					root.node = root.node, model = model, criterion = criterion, shiftCut = shiftCut);

			# Select model with best score according to the specific criterion employed (default aicc)
				#best <- which.min(unlist(lapply(res, "[[", criterion)));
				#best <- res[[which.min(unlist(lapply(res, "[[", criterion)))]];
				fit <- res[[which.min(unlist(lapply(res, "[[", criterion)))]];
				
			# node gives all nodes. doh...
			# this doesn't seem necessary anymore...
				#node <- best$split.at;
				#node <- tail(best$split.at, 1);
				#cut <- best$cut.at;
				
			# here be the problem. this doesn't seem necessary anymore...
			#	fit <- .update.fit.medusa (node = node, z = z, desc = desc, fit = fit, prefit = prefit, num.tips = num.tips,
			#		root.node = root.node, model = model, criterion = criterion, shiftCut = cut);
				#z <- .split.z.at.node.medusa(node = tail(best$split.at, 1), z = z, desc = desc, shiftCut = tail(best$cut.at, 1))$z;
				#step <- rbind(optModel$step, c("add", tail(best$split.at, 1)));
				
				z.best <- .split.z.at.node.medusa(node = tail(fit$split.at, 1), z = z, desc = desc, shiftCut = tail(fit$cut.at, 1))$z;
				step <- rbind(optModel$step, c("add", tail(fit$split.at, 1)));
				
			# Consider parameter removal
				if (stepBack) {
					#backFit <- .back.step.medusa(currentModel = fit, z = z, step = step, model = model, fixPar = fixPar, criterion = criterion);
					#backFit <- .back.step.medusa(currentModel = best, z = z, step = step, model = model, criterion = criterion);
					backFit <- .back.step.medusa(currentModel = fit, z = z.best, step = step, model = model, criterion = criterion);
					fit <- backFit$fit;
					z.best <- backFit$z;
					step <- backFit$step;
				}
				#fit$z <- z;
				fit$step <- step;
			# Compare last accepted model to current best model
				#if (as.numeric(models[[length(models)]][criterion]) - as.numeric(res[[best]][criterion]) < threshold_N) {
				if (as.numeric(optModel[criterion]) - as.numeric(fit[criterion]) < threshold_N) {
					cat("\nNo significant increase in ", criterion, " score. Disregarding subsequent piecewise models.\n\n", sep="");
					done <- TRUE;
					break;
				} else {
				
					if (stepBack && !is.null(backFit$remove)) {.print.removed.shifts(remove=backFit$remove);}
					optModel <- fit;
					z <- z.best;
					
					#models <- c(models, res[best]);
					#fit <- res[[best]]; # keep track of '$split.at' i.e. nodes already considered
					#z <- zz[[i + 1]] <- .split.z.at.node.medusa(node = node.list[best], z = z, desc = desc, shiftCut = fit$cut.at)$z;
					
					#cat("Step ", i+1, ": best likelihood = ", models[[i+1]]$lnLik, "; AICc = ", models[[i+1]]$aicc, "; break at node ",
					#	models[[i+1]]$split.at[i+1], "; model=", models[[i+1]]$model, "; cut=", models[[i+1]]$cut.at, "\n", sep="");
					cat("Step ", i+1, ": lnLik=", round(fit$lnLik, digits=7),
						"; AICc=", fit$aicc, "; shift at node ", tail(fit$split.at,1), "; model=",
						tail(fit$model,1), "; cut=", tail(fit$cut.at,1), "; # shifts=", length(fit$split.at) - 1, "\n", sep="");
					i <- i + 1;
				}
			}
			optModel$z <- z;
			return(optModel);
		}
		
	# temporarily turning off 'stopOnLimit', although probably should be pitched altogether
		#if (stop == "partitions") optModel <- stopOnLimit(fit = fit) else optModel <- stopOnThreshold(fit = fit);
		if (stop == "threshold") optModel <- stopOnThreshold(fit = fit) else optModel <- stopOnLimit(fit = fit);
		
		#modelSummary <- .summary.modelfit.medusa(models = models, phy = phy, threshold = threshold_N);
		modelSummary <- .summary.modelfit.medusa(optModel);

		#if (verbose) {
		#	cat("Model fit summary:", "\n\n", sep = "");
		#	print(modelSummary);
		#	cat("\n");
		#	# if (threshold_N > 0) {
		#		# cat("\nAIC weights are not reported, as they are meaningless when using a threshold criterion.\n");
		#	# }
		#}

		## SUMMARY - this uses a long lost format
		# just output the summary, instead of the function
		# zSummary <- function (id) {
		zSummary <- function () {
			# model.id <- id;
			# #z <- zz[[model.id]];
			# how is opt.model different from modelSummary?
			# opt.model <- data.frame(cut = modelSummary$cut[1:model.id], split = modelSummary$split[1:model.id], models[[model.id]]$par,
				# lnLp = models[[model.id]]$lnLik.part, stringsAsFactors = FALSE);
			opt.model <- data.frame(cut = optModel$cut.at, split = optModel$split.at, optModel$par, lnLp = optModel$lnLik.part, stringsAsFactors = FALSE);
			# break.pts <- opt.model$split;
			break.pts <- as.numeric(optModel$split.at);
			# cut.at <- opt.model$cut;
			cut.at <- optModel$cut.at;
			rownames(z) <- phy$hash[z[, "dec"]]; # identify identical edges among trees
			
			# temp hack. old format used NA for split.at for root
			root <- min(z[,"anc"]);
			break.pts[which(break.pts == root)] <- NA;
			
			# collect shift times
			internal <- is.na(z[, "n.t"]);
			res <- numeric(nrow(z));
			res[] <- NA;
			check <- !is.na(break.pts);
			cut.at <- cut.at[check];
			break.pts <- break.pts[check];
			if (length(break.pts)) {
				for (i in 1:length(break.pts)) {
					idx <- which(z[, "dec"] == break.pts[i]);
					if (cut.at[i] == "node" && internal[idx]) {
						res[idx] <- z[idx, "t.1"];
					} else if (cut.at[i] == "stem" || !internal[idx]) { # pendant node
						res[idx] <- z[idx, "t.0"];
					} else {
						stop("'cut' should either be stem or node");
					}
				}
			}
			anc <- match(z[, "anc"], z[, "dec"]);
			#z <- cbind(z, shift = match(z[, "dec"], opt.model$split), t.shift = res);
			#z <- cbind(z, r = opt.model[z[, "partition"], "r"], epsilon = opt.model[z[, "partition"], "epsilon"]);
			z <- cbind(z, shift = match(z[, "dec"], optModel$split), t.shift = res);
			z <- cbind(z, r = optModel$par[z[, "partition"],"r"], epsilon = optModel$par[z[, "partition"],"epsilon"]);
			z <- cbind(z, ancestral.r = z[anc, "r"], ancestral.epsilon = z[anc, "epsilon"]);
			return(z);
		}

		#model.id <- nrow(modelSummary);
		# z.summary <- zSummary(id = model.id);

		control <- list(stop = stop, threshold = structure(threshold_N, names = criterion), partitions = npartitions);
		#results <- list(control = control, cache = list(desc = desc, phy = phy, richness = richness), model = optModel,
		#	summary = modelSummary, FUN = zSummary);
		results <- list(control = control, cache = list(desc = desc, phy = phy), model = optModel,
			summary = modelSummary, zSummary = zSummary(), medusaVersion = medusaVersion);
		return(results);
	}
	
	## I don't like this format. Why store N (possibly 10,000) copies of richness? Same goes for stop and npartitions
	# seems like hash is ONLY useful with mulitple trees. wtf?!?
	if ("multiPhylo" %in% class(phy)) { ## deal with multiple trees
		#res <- lapply(phyData, function (x) medusa_runner(x$phy, x$richness));
		res <- lapply(phyData$phy, function (x) medusa_runner(phy = x, richness = phyData$richness));
		names(res) <- names(phy);
		class(res) <- c("multimedusa", class(res));
	} else {
		# phyData$phy <- hashes.phylo(phyData$phy);
		res <- medusa_runner(phyData$phy, phyData$richness);
		
		# add profile likelihoods on parameter values
		res$model$prof.par <- .get.profile.likelihoods(res);
		res$summary <- cbind(res$summary, res$model$prof.par);
		
		res$richness <- phyData$richness;
		class(res) <- c("medusa", class(res));
	}

	invisible(res);
}

##############################################################################
##############################################################################


# make sure things are in the correct order and of correct format
.check.richness <- function (phy, richness = NULL) {
	if (is.null(richness)) {
		if ("multiPhylo" %in% class(phy)) {
			phy <- phy[[1]];
		}
		richness <- data.frame(taxon = phy$tip.label, n.taxa = 1);
	} else {
		richness = data.frame(richness, stringsAsFactors = FALSE);
		if (length(richness[1, ]) == 2) {
			if (colnames(richness)[1] != "taxon" || colnames(richness)[2] != "n.taxa") {
				if (class(richness[, 1]) == "factor" & class(richness[, 2]) == "integer") {
					colnames(richness) = c("taxon", "n.taxa");
				} else if (class(richness[, 1]) == "integer" & class(richness[, 2]) == "factor") {
					colnames(richness) = c("n.taxa", "taxon");
				} else {
					stop("'richness' data appear incorrectly formated: see medusa()");
				}
			}
		}
	}
	return(richness);
}

## Function to prune tree using 'richness' information, assumed to have minimally two columns, "taxon" and "n.taxa"
##   Perhaps relax on these column names, may cause too many problems
## May also include 'exemplar' column; in that case, rename relevant tip.label before pruning. this is currently broke.
#prune.tree.merge.data
.treedata.medusa <- function (phy, richness = NULL, ...) {
		
	## MODIFIED -- jme
	if ("multiPhylo" %in% class(phy)) { ## deal with multiple trees
		res = lapply(phy, .treedata.medusa, richness, ...);
		tips = c();
		for (i in 1:length(res)) {
			tips = union(tips, res[[i]]$phy$tip.label);
		}
		nn = sapply(res, function (x) Ntip(x$phy));
		if (!all(nn == length(tips))) {
			stop("'phy' appears to have trees with incompletely overlapping sets of tips");
		}
		trees = lapply(res, "[[", "phy");
		class(trees) = "multiPhylo";
		trees = hashes.phylo(phy = trees, tips = tips) ## FROM PHYLO package [returns trees with $hash object -- a 'key' for each split]
		#for (i in 1:length(res)) res[[i]]$phy = trees[[i]];
		#return(res);
		return(list(phy = trees, richness = richness));
	}
	## END -- jme
	
	# Rename exemplar taxa with taxon name in richness file. Doesn't work in current implementation
	if (!is.null(richness$exemplar)) {
		# Change relevant tip.labels in phy; individual 'exemplar' may be NA, use original tip.label.
		# Ordering in richness file should NOT be assumed to match order of tip.labels
		i.na <- is.na(richness$exemplar);
		phy$tip.label[match(richness$exemplar[!i.na], phy$tip.label)] <- as.character(richness$taxon[!i.na]);
	}

	# checking for typo; if same size, nothing should be dropped
	check <- FALSE;
	if (length(phy$tip.label) == length(richness[, 1])) {
		check <- TRUE;
	}

	# Prune tree down to lineages with assigned richnesses
	temp <- richness[, "n.taxa"];
	names(temp) <- richness[, "taxon"];
	pruned <- treedata(phy, temp, ...); # geiger function calling ape (namecheck)
	if (check) {
		if (length(phy$tip.label) != length(pruned$phy$tip.label)) {
			stop("'richness' data and tip labels in 'phy' do not appear to match fully");
		}
	}
	phy <- pruned$phy;
	rr <- pruned$data;
	richness <- structure(data.frame(taxon = rownames(rr), n.taxa = rr[, 1], row.names = NULL));
	# Check the tree
	# plotNN(phy); # Node numbers (ape-style) plotted

	return(list(phy = phy, richness = richness));
}

## Original default was to fit 20 models (or less if the tree was small).
## Changing to a stop-criterion (stop="model.limit") e.g. when k = n-1 (i.e. when denominator of aicc correction is undefined).
## k <- (3*i-1) # when both birth and death are estimated, where i is the number of piecewise models
## This occurs when i = n/3
## If Yule, max i = n/2
## n <- (2*num.taxa - 1) == (2*length(richness[,1]) - 1) # i.e. total number of nodes in tree (internal + pendant)
## Alternatively use aicc threshold itself as a stopping criterion (stop="threshold").
# AICc = AIC + 2*k*(k+1)/(n-k-1);
.get.max.model.limit <- function (richness, stop, npartitions, model, verbose) {
	samp.size <- (2 * nrow(richness) - 1);
	if (model == "bd" || model == "mixed") {
		max.model.limit <- as.integer(samp.size/3) - ((!(samp.size%%3)) * 1);
	} else {
		max.model.limit <- as.integer(samp.size/2) - ((!(samp.size%%2)) * 1);
	}

	if (stop == "partitions") {
		if (npartitions > max.model.limit) {
			model.limit <- max.model.limit;
			warning("Supplied 'partitions' is in excess of the maximal number that can be considered");
		} else {
			model.limit = npartitions;
		}
	} else {
		model.limit <- max.model.limit;
	}

	if (verbose) {
		cat("Limiting consideration to a maximum of ", model.limit, " piecewise", sep = "");
		if (model == "bd") {
			cat(" birth-death models");
		} else if (model == "mixed") {
			cat(" mixed models");
		} else {
			cat(" pure-birth (Yule) models");
		}
		if (stop == "threshold") {
			cat(" (or until threshold is not satisfied)");
		}
		cat(".\n\n");
	}
	return(model.limit);
}

## Fitted curve from random b-d simulations
## Value corresponds to 95th percentile of AICc(split) - AICc(no-split) for no-split simulations
## x-shifted power function
#get.threshold
#.threshold.medusa <- function (N) {
.threshold.medusa <- function (phy) {
	if ("multiPhylo" %in% class(phy)) {
		phy <- phy[[1]];
	}
	N <- Ntip(phy);
	a <- -35.9410523803326;
	b <- 6.7372587299747;
	c <- -0.100615083407549;
	Offset <- 27.5166786643334;
	y <- a * (N - b)^c + Offset;
	if (y < 0) {
		y <- 0;
	}
	cat("Appropriate AICc threshold for tree of ", N, " tips is: ", y, ".\n\n", sep="");
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
.make.cache.medusa <- function (phy, richness, fx, shiftCut) {
	n.tips <- length(phy$tip.label);
	n.int <- nrow(phy$edge) - n.tips;

	## Ape numbers the tips first
	i.int <- seq_len(n.int);
	interior <- phy$edge[, 2] %in% phy$edge[, 1];
	bt <- branching.times(phy);

	# Consider only internal edges first. may be zero if only 2 tips.
	if (n.int > 0) {
		edges.int <- matrix(phy$edge[interior,], nrow = n.int, ncol = 2);
		colnames(edges.int) <- c("anc", "dec");
		
		t.0 <- bt[match(edges.int[,1], (n.tips + 1):max(edges.int))];
		t.1 <- c(t.0[i.int] - phy$edge.length[interior]);
		
		z.internal <- cbind(edges.int, t.0, t.1, t.len = t.0 - t.1,
			n.0 = rep(1, n.int), n.t = rep(NA, n.int));
	}
	
	# Now, pendant edges;
	edges.pendant <- phy$edge[match(seq_len(n.tips), phy$edge[,2]),];
	colnames(edges.pendant) <- c("anc", "dec");

	t.0 <- bt[match(edges.pendant[, 1], (n.tips + 1):max(edges.pendant))];
	t.1 <- rep(0, n.tips);
	# cannot assume richness ordering necessarily matches that of tip labels
	ext.richness <- richness$n.taxa[match(phy$tip.label, richness$taxon)];

	z.pendant <- cbind(edges.pendant, t.0, t.1, t.len = t.0 - t.1,
		n.0 = rep(1, n.tips), n.t = ext.richness);
	
	if (n.int > 0) {
		z <- rbind(z.internal, z.pendant);
	} else { # case with only 2 pendant edges.
		z <- z.pendant;
	}
	
	z <- cbind(z, partition = rep(1, length(z[, 1]))) # Stores piecewise model structure
	rownames(z) <- NULL;

	# Used for identifying descendant nodes below i.e. tracking breakpoints
	all.edges <- as.matrix(z[, c("anc", "dec")]);
	desc.stem <- list();
	desc.node <- list();
	
	# only calculate what is needed given value of shiftCut
	if (shiftCut == "both" || shiftCut == "stem") {
		desc.stem <- fx(seq_len(max(all.edges)), .descendants.cutAtStem.idx, all.edges = all.edges);
	}
	if (shiftCut == "both" || shiftCut == "node") {
		if (!is.null(desc.stem)) {
			root <- min(z[,"anc"]);
			desc.node <- fx(desc.stem, .strip.stem);
			desc.node[root] <- desc.stem[root];
		} else {
			desc.node <- fx(seq_len(max(all.edges)), .descendants.cutAtNode.idx, all.edges = all.edges);
		}
	}
	
	# ?!? wtf is this? Get rid of it?
	for (i in 1:length(desc.node)) {
		if (length(desc.node[[i]]) == 0) { # tips
			desc.node[[i]] <- .descendants.cutAtStem.idx(node.list = i, all.edges = all.edges);
		}
	}

	tips <- .cache.descendants(phy)$tips;
	res <- list(z = z, desc.stem = desc.stem, desc.node = desc.node, tips = tips);
	return(res);
}

## This generates the indices of all descendants of a node, using ape's edge matrix.
## Deals with row numbers of the edge matrix rather than node numbers of the tree.
.descendants.cutAtStem <- function (node, all.edges) {
	ans <- numeric();
	ans <- node;
	repeat {
		node <- all.edges[all.edges[, 1] %in% node, 2];
		if (length(node) > 0) {
			ans <- c(ans, node);
		} else {
			break;
		}
	}
	return(unlist(ans));
}

## The function 'descendants' returns the indices of all descendants within the edge matrix.
.descendants.cutAtStem.idx <- function (node.list, all.edges) {
	which(all.edges[, 1] == node.list | all.edges[, 2] %in% .descendants.cutAtStem(node.list, all.edges));
}

## This generates the indices of all descendants of a node, using ape's edge matrix.
## Deals with row numbers of the edge matrix rather than node numbers of the tree.
.descendants.cutAtNode <- function (node, all.edges) {
	ans <- numeric();
	repeat {
		node <- all.edges[all.edges[, 1] %in% node, 2]
		if (length(node) > 0) {
			ans <- c(ans, node);
		} else {
			break;
		}
	}
	return(unlist(ans));
}

## The function 'descendants' returns the indices of all descendants within the edge matrix.
.descendants.cutAtNode.idx <- function (node.list, all.edges) {
	which(all.edges[, 1] == node.list | all.edges[, 2] %in% .descendants.cutAtNode(node.list, all.edges));
}

# Remove stem node from previously calculated set of stem descendants; about a billion times faster than descendants.cutAtNode
.strip.stem <- function (x) {
	y <- unlist(x);
	return(y[-1]);
}


## Needed for determining whether nodes are virgin nodes
# ?!? this isn't even used anywhere...
#get.num.tips
.ntips <- function (node, phy) {
	n <- length(tips(phy, node));
	return(n);
}

## Only used for base model
#.fit.base.medusa <- function (z, sp, model, fixPar, criterion) {
.fit.base.medusa <- function (z, sp, model, criterion) {
	#fit <- .get.optimal.model.flavour(z = z, sp = sp, model = model, fixPar = fixPar, criterion = criterion);
	fit <- .get.optimal.model.flavour(z = z, sp = sp, model = model, criterion = criterion);
	model.fit <- .calculate.modelfit.medusa(fit = fit, z = z);
	steps <- c("add", as.numeric(min(z[,"anc"])));
	
	#return(list(par=matrix(fit$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik.part=fit$lnLik, lnLik=fit$lnLik,
	#	split.at=min(z[,"anc"]), aic=round(model.fit$aic, digits=7), aicc=round(model.fit$aicc, digits=7), num.par=model.fit$k,
	#	cut.at="node", model=fit$model, z=z, step=matrix(steps, nrow=1, dimnames=list(NULL,c("step", "node")))));
	return(list(par=matrix(fit$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik.part=fit$lnLik, lnLik=fit$lnLik,
		split.at=min(z[,"anc"]), aic=round(model.fit$aic, digits=7), aicc=round(model.fit$aicc, digits=7), num.par=model.fit$k,
		cut.at="node", model=fit$model, step=matrix(steps, nrow=1, dimnames=list(NULL,c("step", "node")))));
}

## When model == mixed, fit both and find optimal flavour
#.get.optimal.model.flavour <- function (z, sp, model, fixPar, criterion) {
.get.optimal.model.flavour <- function (z, sp, model, criterion) {
	fit.bd <- NULL;
	fit.yule <- NULL;
	fit <- NULL;
	
	if (model == "yule" | model == "mixed") {
		fit.yule <- .fit.partition.medusa(z = z, sp = sp, model = "yule");
		fit.yule$model <- "yule";
	}
	if (model == "bd" | model == "mixed") {
		if (is.na(sp[2])) {sp[2] <- 0.5;}
		fit.bd <- .fit.partition.medusa(z = z, sp = sp, model = "bd");
		fit.bd$model <- "bd";
	}
	# if (model != "mixed" && model != "bd" && model != "yule") { # i.e. the constrained models
		# fit <- medusaMLFitPartition(z=z, sp=sp, model=model, fixPar=fixPar);
		# fit$model <- model;
		# return(fit);
	# }
## Figure out which model fits best
	if (is.null(fit.bd)) {
		fit <- fit.yule;
	} else if (is.null(fit.yule)) {
		fit <- fit.bd;
	} else {
## Considering both models
		fit <- .get.best.partial.model(fit1 = fit.yule, fit2 = fit.bd, z = z, criterion = criterion);
	}
	return(fit);
}

## Used for comparing models fit to the same partition
.get.best.partial.model <- function (fit1, fit2, z, criterion) {
	n <- (length(z[,1]) + 1); # the number of nodes involved
## Add '1' to parameters to account for break
	k1 <- 1 + sum(!is.na(fit1$par));
	k2 <- 1 + sum(!is.na(fit2$par));
	
	if (n - k1 <= 1 || n - k2 <= 1) { # deals with single edges, where AICc correction becomes undefined. use AIC.
		if (.get.AIC(fit1$lnLik, k1) < .get.AIC(fit2$lnLik, k2)) {
			return(fit1);
		} else {
			return(fit2);
		}
	} else {
		if (criterion == "aicc") {
			if (.get.AICc(fit1$lnLik, k1, n) < .get.AICc(fit2$lnLik, k2 ,n)) {
				return(fit1);
			} else {
				return(fit2);
			}
		} else {
			if (.get.AIC(fit1$lnLik, k1) < .get.AIC(fit2$lnLik, k2)) {
				return(fit1);
			} else {
				return(fit2);
			}
		}
	}
}





# replace these at some point with a general function
.get.AIC <- function (lnLik, k) {
	return(-2 * lnLik + 2*k);
}

.get.AICc <- function (lnLik, k, n) {
	return(-2 * lnLik + 2*k*n/(n-k-1));
}



#.prefit.tip.medusa <- function (node, z, sp, model, fixPar, criterion) {
.prefit.tip.medusa <- function (node, z, sp, model, criterion) {
	z.tip <- z[z[,"dec"] == node,,drop=FALSE];
	
	if (all(z.tip[,"n.t"] == 1) && (model == "yule" || model == "mixed")) {
		return(list(par=c(0, NA), lnLik=0, model="yule"));
	}
	# tips are always better fit by yule. don't bother with BD.
	#fit <- .get.optimal.model.flavour(z=z.tip, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
	if (model == "yule" || model == "mixed") {
		return(.get.optimal.model.flavour(z=z.tip, sp=sp, model="yule", criterion=criterion));
	} else { # at the moment, only BD will get through. but eventually constrained models also
		return(.get.optimal.model.flavour(z=z.tip, sp=sp, model=model, criterion=criterion));
	}
}



	
.prefit.medusa <- function (node, z, desc, sp, model, shiftCut, criterion) {
	zz <- .split.z.at.node.medusa(node = node, z = z, desc = desc, shiftCut = shiftCut, extract = TRUE)$z;
	#fit <- .get.optimal.model.flavour(z=z.node, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
	fit <- .get.optimal.model.flavour(z = zz, sp = sp, model = model, criterion = criterion);
	return(fit);
}



## Pre-fit values for pendant edges; DON'T recalculate later; should account for ~25% of all calculations
## Also cache values for virgin nodes; useful until subsetted.
## shiftCut can only be "stem" or "node" (not "both"), as both are evaluated separately
#medusa.ml.prefit
# .prefit.medusa <- function (node, z, desc, sp, model, shiftCut, criterion) {
	# fitted.bd <- NULL;
	# fitted.yule <- NULL;
	# #   fit <- NULL;
	
	# ## if model is mixed, grab the optimally fitted one, drop the other; it is cutAtStem that matters, as it is the sum of 2 break likelihoods
	# ## optimal alone may be different than optimal in tandem
	
	
	# ## update to use medusaFitOptimal
	# if (shiftCut == "stem") {
		# if (model == "bd" || model == "mixed") {
			# obj <- .split.z.at.node.medusa(node = node, z = z, desc = desc, shiftCut = "stem");
			# z.bd.stem <- obj$z;
			# fitted.bd <- .fit.partition.medusa(partition = 2, z = z.bd.stem, sp = sp, model = "bd");
			# fitted.bd$model <- "bd";
			# fitted.bd$cut.at <- "stem";
		# }
		# if (model == "yule" || model == "mixed") {
			# obj <- .split.z.at.node.medusa(node = node, z = z, desc = desc, shiftCut = "stem");
			# z.yule.stem <- obj$z;
			# fitted.yule <- .fit.partition.medusa(partition = 2, z = z.yule.stem, sp = sp, model = "yule");
			# fitted.yule$model <- "yule";
			# fitted.yule$cut.at <- "stem";
		# }
	# } else if (shiftCut == "node") {
		# if (model == "bd" || model == "mixed") {
			# obj <- .split.z.at.node.medusa(node = node, z = z, desc = desc, shiftCut = "node");
			# z.bd.node <- obj$z;
			# fitted.bd <- .fit.partition.medusa(partition = 2, z = z.bd.node, sp = sp, model = "bd");
			# fitted.bd$model <- "bd";
			# fitted.bd$cut.at <- "node";
		# }
		# if (model == "yule" || model == "mixed") {
			# obj <- .split.z.at.node.medusa(node = node, z = z, desc = desc, shiftCut = "node");
			# z.yule.node <- obj$z;
			# fitted.yule <- .fit.partition.medusa(partition = 2, z = z.yule.node, sp = sp, model = "yule");
			# fitted.yule$model <- "yule";
			# fitted.yule$cut.at <- "node";
		# }
	# }
	# ## Check which flavour of model fits best
	# if (is.null(fitted.bd) & !is.null(fitted.yule)) {
		# return(fitted.yule);
	# } else if (is.null(fitted.yule) & !is.null(fitted.bd)) {
		# return(fitted.bd);
	# } else {
		# ## Dealing with a 'mixed' model here; need to consider number of parameters.
		# bd.model.fit <- .calculate.modelfit.medusa(fit = fitted.bd, z = z)
		# yule.model.fit <- .calculate.modelfit.medusa(fit = fitted.yule, z = z)
		# element <- NULL
		# if (criterion == "aic") {
			# element <- 1;
		# } else {
			# element <- 2;
		# }
		# if (is.nan(yule.model.fit[[element]])) {
			# return(fitted.bd);
		# } else if (bd.model.fit[[element]] < yule.model.fit[[element]]) {
			# return(fitted.bd);
		# } else {
			# return(fitted.yule);
		# }
	# }
# }


## sp = initializing values for r & epsilon
## Default values should never be used (except for first model), as the values from the previous model are passed in
#medusa.ml.fit.partition
#.fit.partition.medusa <- function (partition, z, sp = c(0.1, 0.05), model) {
.fit.partition.medusa <- function (z, sp = c(0.1, 0.05), model) {
	# Construct likelihood function:
	#lik <- .lik.partition.medusa(partition = (z[z[, "partition"] == partition, , drop = FALSE]), model = model);
	lik <- .lik.partition.medusa(partition = z, model = model);
	foo <- function (x) {
		-lik(pars = exp(x));
	} # work with parameters in log-space to preserve precision

	if (model == "bd") {
		fit <- optim(fn = foo, par = log(sp), method = "N"); # last argument connotes maximization
		return(list(par = exp(fit$par), lnLik = -fit$value));
	} else {
		fit <- optimize(f = foo, interval = c(-25, 1));
		par <- c(exp(fit$minimum), NA);
		return(list(par = par, lnLik = -fit$objective));
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
.split.z.at.node.medusa <- function (node, z, desc, shiftCut, extract=FALSE) {
	descendants <- NULL;
	if (shiftCut == "stem") {
		descendants <- desc$stem;
	} else {
		descendants <- desc$node;
	}
	part <- z[, "partition"];
	base <- min(part[z[, 1] == node | z[, 2] == node]);
	tag <- max(part) + 1;
	i <- descendants[[node]];
	idx <- i[part[i] == base];
	z[idx, "partition"] <- tag;
	#z[which(z["dec"] == node), "partition"] <- tag; # Possible to have several edges to consider
	
	if (extract) {z <- z[idx,,drop = FALSE];}
	
	return(list(z = z, affected = c(unique(part[idx]), tag)));
}

## 'fit' contains parameter values from previous model, used to initialize subsequent model.
## Pass in pre-fitted values for pendant edges and virgin nodes (in 'prefit'); DON'T recalculate.
## Need to consider the possibility of birth-death, yule, or mixed models.
## Need to consider where shft is placed (shiftCut). Placement affects not only new clade, but
## also the size of the split clade. Only relevant if shiftCut = "both".
## fit1 model is already logged; only need to record fit2 model, and only non-prefitted nodes
#medusa.ml.update
.update.fit.medusa <- function (node, z, desc, fit, prefit, num.tips, root.node, model, criterion, shiftCut) {
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

	if (shiftCut == "stem" || shiftCut == "both") {
	## First, diminshed clade
		obj.stem <- .split.z.at.node.medusa(node = node, z = z, desc = desc, shiftCut = "stem");
		z.stem <- obj.stem$z;
		aff <- obj.stem$affected;
		
	## Ensure that neither partition is empty; can occur with "node" or "both" cutting. If so, kill it.
		if (sum(z.stem[,"partition"] == aff[1]) == 0 || sum(z.stem[,"partition"] == aff[2]) == 0) {
			fit$lnLik <- -Inf;
			fit$aic <- Inf;
			fit$aicc <- Inf;
			return(fit);
		}
		
	## Everything is cool; proceed.
		sp <- op[aff[1], ]; # Use previously fit parameter values from clade that is currently being split
		
		dimClade <- z.stem[z.stem[,"partition"] == aff[1],,drop = FALSE];
		
	## first, consider diminshed clade. may result in a clade that has been cached previously
	## check if diminished clade conforms to a cached clade
		x <- which(dimClade[,"anc"] == min(dimClade[,"anc"]));
		if (length(x) == 1) { # a cut-at-stem scenario
			y <- as.numeric(dimClade[x, "dec"]);
			if (length(unique(dimClade[(dimClade[,"dec"] < root.node),"dec"])) == num.tips[[y]]) {
				if (y < root.node) {
					fit1.stem <- prefit$tips[[y]];
				} else {
					fit1.stem <- prefit$virgin.nodes$stem[[y - root.node]];
				}
			}
		}
		
		if (is.null(fit1.stem)) { # not cached
			fit1.stem <- .get.optimal.model.flavour(z = dimClade, sp = sp, model = model, criterion = criterion);
			#if (model == "mixed") { ## In mixed models, want to conserve flavour of previously fit model (right?). No!
			#	if (sum(!is.na(sp)) < 2) { # yule
			#		fit1.stem <- .fit.partition.medusa(partition = aff[1], z = z.stem, sp = sp, model = "yule");
			#	} else {
			#		fit1.stem <- .fit.partition.medusa(partition = aff[1], z = z.stem, sp = sp, model = "bd");
			#	}
			#} else {
			#	fit1.stem <- .fit.partition.medusa(partition = aff[1], z = z.stem, sp = sp, model = model);
			#	fit1.stem$model <- model;
			#}
		}
		
	## Second, new clade
		if (node < root.node) { # tip, already calculated
			fit2.stem <- prefit$tips[[node]];
		} else if (length(unique(z.stem[(z.stem[, "partition"] == aff[2] & z.stem[, "dec"] < root.node), "dec"])) == num.tips[[node]]) { # virgin node, already calculated
			fit2.stem <- prefit$virgin.nodes$stem[[node - root.node]];
		} else { # novel shift
			newClade <- z.stem[z.stem[,"partition"] == aff[2],,drop = FALSE]; # only subset if necessary
			fit2.stem <- .get.optimal.model.flavour(z=newClade, sp=sp, model=model, criterion=criterion);
			
			# fit2.stem.bd <- NULL;
			# fit2.stem.yule <- NULL;

			# if (model == "yule" || model == "mixed") {
				# #fit2.stem.yule <- .fit.partition.medusa(aff[2], z.stem, sp = sp, model = "yule");
				# fit2.stem.yule <- .fit.partition.medusa(newClade, sp = sp, model = "yule");
				# fit2.stem.yule$model <- "yule";
			# }
			# if (model == "bd" || model == "mixed") {
				# if (is.na(sp[2])) {
					# sp[2] <- 0.5;
				# }
				# #fit2.stem.bd <- .fit.partition.medusa(aff[2], z.stem, sp = sp, model = "bd");
				# fit2.stem.bd <- .fit.partition.medusa(newClade, sp = sp, model = "bd");
				# fit2.stem.bd$model <- "bd";
			# }
			# ## Figure out which model fits best
			# if (is.null(fit2.stem.bd)) {
				# fit2.stem <- fit2.stem.yule;
			# } else if (is.null(fit2.stem.yule)) {
				# fit2.stem <- fit2.stem.bd;
			# } else {
				# ## Considering both places for a shift
				# fit2.stem.bd.val <- .calculate.modelfit.medusa(fit = fit2.stem.bd, z = z);
				# fit2.stem.yule.val <- .calculate.modelfit.medusa(fit = fit2.stem.yule, z = z);

				# if (criterion == "aic") {
					# element <- 1;
				# } else {
					# element <- 2;
				# }

				# if (fit2.stem.bd.val[[element]] < fit2.stem.yule.val[[element]]) {
					# fit2.stem <- fit2.stem.bd;
				# } else {
					# fit2.stem <- fit2.stem.yule;
				# }
			# }
		}
	}
	
	if ((shiftCut == "node" || shiftCut == "both")  && (node > root.node)) {
		## First, diminshed clade
		obj.node <- .split.z.at.node.medusa(node = node, z = z, desc = desc, shiftCut = "node");
		z.node <- obj.node$z;
		aff <- obj.node$affected;
		sp <- op[aff[1], ]; # Use previously fit parameter values from clade that is currently being split
		
		dimClade <- z.node[z.node[,"partition"] == aff[1],,drop = FALSE];
		
		fit1.node <- .get.optimal.model.flavour(z = dimClade, sp = sp, model = model, criterion = criterion);
		#if (model == "mixed") { ## In mixed models, want to conserve flavour of previously fit model (right?) No!
		#	if (sum(!is.na(sp)) < 2) { # yule
		#		fit1.node <- .fit.partition.medusa(partition = aff[1], z = z.node, sp = sp, model = "yule");
		#		fit1.node$model <- "yule";
		#	} else {
		#		fit1.node <- .fit.partition.medusa(partition = aff[1], z = z.node, sp = sp, model = "bd");
		#		fit1.node$model <- "bd";
		#	}
		#} else {
		#	fit1.node <- .fit.partition.medusa(partition = aff[1], z = z.node, sp = sp, model = model);
		#	fit1.node$model <- model;
		#}
		
		## Second, new clade
		#if (node < root.node) { # not possible
		#	fit2.node <- prefit$tips[[node]];
		#}
		if (length(unique(z.node[(z.node[, "partition"] == aff[2] & z.node[, "dec"] < root.node), "dec"])) == num.tips[[node]]) { # virgin node
			fit2.node <- prefit$virgin.nodes$node[[node - root.node]];
		} else { # novel shift
			newClade <- z.node[z.node[,"partition"] == aff[2],,drop = FALSE];
			fit2.node <- .get.optimal.model.flavour(z=newClade, sp=sp, model=model, criterion=criterion);
			# fit2.node.bd <- NULL;
			# fit2.node.yule <- NULL;

			# if (model == "yule" || model == "mixed") {
				# fit2.node.yule <- .fit.partition.medusa(aff[2], z.node, sp, model = "yule");
				# fit2.node.yule$model <- "yule";
			# }
			# if (model == "bd" || model == "mixed") {
				# if (is.na(sp[2])) {
					# sp[2] <- 0.5;
				# }
				# fit2.node.bd <- .fit.partition.medusa(aff[2], z.node, sp, model = "bd");
				# fit2.node.bd$model <- "bd";
			# }
			# ## Figure out which model fits best
			# if (is.null(fit2.node.bd)) {
				# fit2.node <- fit2.node.yule;
			# } else if (is.null(fit2.node.yule)) {
				# fit2.node <- fit2.node.bd;
			# } else {
				# ## Considering both places for a shift
				# fit2.node.bd.val <- .calculate.modelfit.medusa(fit = fit2.node.bd, z = z);
				# fit2.node.yule.val <- .calculate.modelfit.medusa(fit = fit2.node.yule, z = z);

				# if (criterion == "aic") {
					# element <- 1;
				# } else {
					# element <- 2;
				# }

				# if (fit2.node.bd.val[[element]] < fit2.node.yule.val[[element]]) {
					# fit2.node <- fit2.node.bd;
				# } else {
					# fit2.node <- fit2.node.yule;
				# }
			# }
		}
	}

	## Now, figure out which shift position is optimal
	if (is.null(fit2.node)) {
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
		stem.par <- rbind(fit1.stem$par, fit2.stem$par);
		stem.val <- list(lnLik = stem.lik, par = stem.par);
		stem.fit <- .calculate.modelfit.medusa(fit = stem.val, z = z);

		node.lik <- (fit1.node$lnLik + fit2.node$lnLik);
		node.par <- rbind(fit1.node$par, fit2.node$par);
		node.val <- list(lnLik = node.lik, par = node.par);
		node.fit <- .calculate.modelfit.medusa(fit = node.val, z = z);

		# element <- NULL;
		# if (criterion == "aic") {element <- 1;} else {element <- 2;}

		if (stem.fit[[criterion]] < node.fit[[criterion]]) {
			fit1 <- fit1.stem;
			fit2 <- fit2.stem;
			cut.at <- "stem";
		} else {
			fit1 <- fit1.node;
			fit2 <- fit2.node;
			cut.at <- "node";
		}
	}
	op[aff[1], ] <- fit1$par # Replace parameters with new values for diminished clade
	fit$model[aff[1]] <- fit1$model; # update altered model

	fit$par <- rbind(op, fit2$par);
	fit$lnLik.part[aff] <- c(fit1$lnLik, fit2$lnLik); # Replace parameters with new values for diminished clade
	fit$split.at <- c(fit$split.at, node);
	fit$lnLik <- sum(fit$lnLik.part);

	model.fit <- .calculate.modelfit.medusa(fit = fit, z = z);

	fit$aic <- model.fit$aic;
	fit$aicc <- model.fit$aicc;
	fit$num.par <- model.fit$k;
	#fit$cut.at <- cut.at;
	#fit$model <- fit2$model;
	
	fit$cut.at <- c(fit$cut.at, cut.at);
	fit$model <- c(fit$model, fit2$model);
	
	return(fit);
}






### NEED TO PUT IN ALL OTHER MODELS HERE!!! ###

## make.lik.medusa.part: generate a likelihood function for a single partition.
#make.lik.medusa.part
## lots of alternative constrained models are missing - add back in
.lik.partition.medusa <- function (partition, model) {
	# Handle internal and pendant edges separately
	is.int <- is.na(partition[, "n.t"]);
	is.pend <- !is.int;

	n.int <- sum(is.int);
	n.pend <- sum(is.pend);

	if (n.int + n.pend != length(partition[, 1])) {
		stop("You messed up, yo.");
	}

	## Internal and pendant calculations differ; split'em up
	int <- partition[is.int, , drop = FALSE];
	pend <- partition[is.pend, , drop = FALSE];

	sum.int.t.len <- sum(int[, "t.len"]); # Simply sum all internal edges
	int.t.0 <- int[, "t.0"];

	# 'n.0' = Foote's 'a', initial diversity; 'n.t' = Foote's 'n', final diversity
	pend.n.0 <- pend[, "n.0"]; # Foote's 'a': initial diversity
	pend.n.t <- pend[, "n.t"]; # Foote's 'n': final diversity
	pend.t.len <- pend[, "t.len"];

	# User may pass in epsilon; don't change it, just estimate r
	f <- function (pars) {
		if (model == "bd") {
			r <- pars[1];
			epsilon <- pars[2];

			if (r < 0 || epsilon <= 0 || epsilon >= 1) {
				return(-Inf);
			}
		} else {
			r <- pars[1];
			epsilon <- 0;
			if (r < 0) {
				return(-Inf);
			}
		}
		l.int <- numeric();
		l.pend <- numeric();

		if (n.int == 0) {
			l.int <- 0;
		} else {
			## Likelihood of internal edges from Rabosky et al. (2007) equation (2.3):
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
		}

		if (n.pend == 0) {
			l.pend <- 0;
		} else {
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
# i.pend.n.t.1 <- which(pend.n.t == 1)   # calculations even simpler: log(1-B)
# i.pend.n.t.n1 <- which(pend.n.t != 1)

			ert <- exp(r * pend.t.len);
			B <- (ert - 1)/(ert - epsilon); # Equivalently: B <- (bert - b) / (bert - d)

			l.pend <- sum(log(1 - B) + (pend.n.t - 1) * log(B));
		}
		return(l.int + l.pend);
	}
}


## 'fit' contains '$par' and '$lnlik'
#calculate.model.fit
.calculate.modelfit.medusa <- function (fit, z) {
	## Sample size taken (for now) as the total num.nodes in the tree (internal + pendant)
	# num.nodes = (2*length(phy$tip.label) - 1) == (2*length(richness[,1]) - 1) == length(z[,1]) + 1
	#   n <- (length(z[,1]) + 1) + sum(!is.na(z[,"n.f"]));

	# Since each edge defines a node (i.e. an 'observation'), need only add root node as final obervation
	n <- (length(z[, 1]) + 1);

	# Includes both formal parameters AND number of breaks. Note: first model does not involve a break.
	## Models where all parameters are estimated (i.e. BD model):
	# 2 parameters for base model (no breakpoint) + 3 parameters (r, eps, breakpoint) for each subsequent model

	# Determine number of piecewise models currently involved
	if (length(fit$par) < 3) { # i.e. base model
		num.models <- 1;
	} else {
		num.models <- length(fit$par[, 1]);
	}

	# Updated for more general models: check how many parameter values != NA
	#   k <- 2 + (3 * (num.models - 1))
	k <- sum(!is.na(fit$par)) + (num.models - 1); # number of estimated parameters + number of breaks

	ll <- list(k = k, lnL = fit$lnLik);
	return(.aic(ll, n));
	#return(.aic(ll$lnL, n, k));
	#   lnLik <- fit$lnLik;
	#   aic <- (-2 * lnLik) + (2*k);
	#   aicc <- aic + 2*k*(k+1)/(n-k-1);
	#   model.fit <- c(aic, aicc, k);
	#   return(model.fit);
}

## Prints out a table of likelihoods, and parameters.
## because a threshold is used, aic weights are meaningless.
.summary.modelfit.medusa <- function (optModel) {
	modelSize <- length(optModel$split.at);
	
	summ <- data.frame(cbind(seq(1:modelSize), optModel$split.at, optModel$cut.at, optModel$model,
		signif(optModel$lnLik.part, digits=7)), signif(optModel$par, digits=6), stringsAsFactors = FALSE);
	colnames(summ) <- c("Model.ID", "Shift.Node", "Cut.At", "Model", "Ln.Lik.part", "r", "epsilon");
	
	return(summ);
}




# WTF? This is ancient code...
## Prints out a table of likelihoods, parameters, and aic scores
# calculate.model.fit.summary
# .summary.modelfit.medusa <- function (models, phy, threshold = 0) {
	# tmp <- matrix(nrow = (length(models)), ncol = 6);
	# colnames(tmp) <- c("partitions", "split", "lnL", "k", "aic", "aicc");
	# # c("N.Models", "Break.Node", "Ln.Lik", "N.Param", "aic", "aicc")
	
	# w.aic <- numeric(length(models));
	# w.aicc <- numeric(length(models));
	# cut.at <- character(length(models));

	# for (i in 1:length(tmp[, 1])) {
		# tmp[i, ] <- c(i, as.integer(models[[i]]$split.at[i]), models[[i]]$lnLik, models[[i]]$num.par, models[[i]]$aic, models[[i]]$aicc);
		# cut.at[i] <- models[[i]]$cut.at;
	# }

	# all.res <- data.frame(tmp);
	# all.res[1, 2] <- NA; # root node for base model
	# all.res$cut <- as.character(cut.at);
	# all.res <- all.res[, c("partitions", "split", "cut", "lnL", "k", "aic", "aicc")];
	
	# # these are meaningless
	# if (threshold == 0) {
		# w.aic <- round(aicw(all.res$aic), digits = 5);
		# w.aicc <- round(aicw(all.res$aicc), digits = 5);
		# all.res$w.aic <- w.aic$w;
		# all.res$w.aicc <- w.aicc$w;
	# }
	# #   if (plot)
	# #   {
	# # 		dev.new();
	# # 		.plot.modelfit.medusa(all.res, ...);
	# # 		if (!is.null(fig.title)) {title(main=fig.title, cex.main=0.75);}
	# #   }
	# return(all.res);
# }


## Self explanatory. Don't use these.
## These are meaningless when using a threshold criterion
#calculate.model.weights
aicw <- function (x) {
	aic <- x;
	best <- min(aic);
	delta <- aic - best;
	sumDelta <- sum(exp(-0.5 * delta));
	w <- (exp(-0.5 * delta)/sumDelta);

	results <- data.frame(fit = aic, delta = delta, w = w);
	rownames(results) <- names(aic);

	return(results);
}

# not useful anymore
#summary.medusaRAW <- function (object, criterion = c("aicc", "aic"), model = NULL, ...) {
.summary.medusa <- function (object, ...) {
	#    function (results, modelNum=NULL, cutoff="threshold", criterion="aicc", plotTree=TRUE, time=TRUE, node.labels=TRUE, cex=0.5, plotSurface=FALSE, n.points=100, ...)
	# Desirables:
#  1. table listing parameter values of selected model <- no longer relevant
#  2. list parameters of base model <- no longer relevant
#  3. tree printed with colour-coded edges, node labels to indicate split position(s)
#  4. plot likelihood surface <- where did this go?!?!?!?!?!?

	# Extract constituent components from results
	#    ct=list(model=NULL, criterion=c(NA, "aic", "aicc"), threshold=0)
#    ct[names(control)]=control

	#criterion <- match.arg(criterion, c("aicc", "aic"));

	results <- object;
	#fit <- results$models;
	fit <- results$model;
	phy <- results$cache$phy;
	desc <- results$cache$desc;
	modelSummary <- results$summary;
	#cutAt <- as.character(modelSummary$cut);
	cutAt <- as.character(modelSummary$Cut.At);

	# fit; phy; z; desc; modelSummary; threshold;
	
	# First, determine which model is desired <- this is deprecated, since only a single (optimal) model is held.
	# model.chooser <- function (criterion = c("aic", "aicc"), threshold = "auto") {
		# if (threshold == "auto") {
			# threshold <- .threshold.medusa(Ntip(phy));
		# } else if (!is.numeric(threshold)) {
			# stop("'threshold' must be numeric");
		# }
		# threshold <- abs(threshold);

		# model.id <- 1;
		# while (1) {
			# if ((model.id + 1) > length(fit)) {
				# break;
			# }
			# if ((fit[[model.id]][[criterion]] - fit[[model.id + 1]][[criterion]]) < threshold) {
				# break;
			# }
			# model.id <- model.id + 1;
		# }
		# structure(threshold, names = model.id);
	# }

	# model.id <- 0
	# if (!is.null(model)) {
		# xx <- list(...);
		# if ("threshold" %in% names(xx)) {
			# warning("'threshold' has been ignored");
		# }
		# model.id <- model;
		# threshold <- 0;
	# } else { # Find best model using some criterion (threshold or user-defined)
		# # if (cutoff != "threshold") {threshold <- cutoff}
		# # else {cat("\nSelecting model based on corrected threshold (improvement in information theoretic score of ", ct$threshold, " units).\n", sep="");}
		# tmp <- model.chooser(criterion = criterion, ...);
		# model.id <- as.integer(names(tmp));
		# threshold <- as.numeric(tmp);
	# }

	#names(threshold) <- criterion;
	#chosen <- ifelse(1:length(fit) %in% model.id, "*", "");
	#modelSummary$chosen <- chosen;
	#z <- object$FUN(model.id);
	#z <- object$FUN();
	z <- object$zSummary;
	attr(z, "phylo") <- phy;
	attr(z, "summary") <- modelSummary;
	#attr(z, "threshold") <- threshold;
	class(z) <- c("medusa", class(z));
	return(z);
}

print.medusa <- function (x, ...) {
	cat("\nOptimal medusa model for tree with ", length(x$cache$phy$tip.label), " taxa.\n\n", sep="");
	print(x$summary);
	cat("\n");
	cat("95% confidence intervals on parameter values calculated from profile likelihoods\n");
}

#plot.medusa <- function (x, partitions = list(cex = 2, bg = "gray", alpha = 0.75, col = "black", lwd = 1), ...) {
plot.medusa <- function (x, cex = 0.5, time = TRUE, bg = "gray", alpha = 0.75, col = "black", lwd = 1, ...) {
	#z <- x;
	z <- x$zSummary;
	#st <- list(cex = 2, bg = "gray", alpha = 0.75, col = "black", lwd = 1);
	#st[names(partitions)] <- partitions;
	#shift <- st;
	shift <- list(cex = cex, bg = bg, alpha = alpha, col = col, lwd = lwd);
	#phy <- attr(z, "phylo");
	phy <- x$cache$phy;
	# par=match.arg(par, c("r", "epsilon"));
	# bb=matrix(z[,par], nrow=1);
	# bb=rbind(bb, bb);
	# colnames(bb)=rownames(z);
	# shifts=max(z[,"partition"]);
	# yy=.branchcol.plot(phy, bb, plot=FALSE, colors=list(branches=min(c(256, 4*shifts+1)), legend=7, missing=1))$col;
	plot(phy, cex=0.5, ...);
	if (time) axisPhylo();
	shifts <- z[(idx <- !is.na(z[, "shift"])), "dec"];
	if (length(shifts)) {
		ss <- z[idx, "shift"];
		ww <- character(nrow(phy$edge));
		ww[match(shifts, phy$edge[, 2])] <- ss;
		xx <- numeric(nrow(phy$edge));
		xx[phy$edge[, 2] %in% shifts] <- 1;
		edgelabels.auteur(NULL, frame = "circle", cex = ifelse(xx == 1, shift$cex, 1e-10), pch = ifelse(xx == 1, 21, NA),
		bg = .transparency(shift$bg, shift$alpha), col = shift$col, lwd = shift$lwd);
		edgelabels.auteur(ww, frame = "none", cex = ifelse(xx == 1, shift$cex/3, 1e-10), pch = NA, col = shift$col);
	}
}

# long deprecated...
# plot.medusaRAW = function (x, col = c(aic = "black", aicc = "blue"), ...) {
	# all.res <- x$summary;
	# ylim <- c(min(all.res[, "aic"], all.res[, "aicc"]), max(all.res[, "aic"], all.res[, "aicc"]));
	# plot(all.res[, "partitions"], all.res[, "aicc"], xlab = "piecewise models", ylab = "model fit", ylim = ylim, type = "l", col = col[["aicc"]], xaxt = "n", ...);
	# axis(1, at = all.res[, "partitions"], labels = all.res[, "partitions"], tcl = NA);
	# points(all.res[, "partitions"], all.res[, "aicc"], col = col[["aicc"]], pch = 21, bg = "white");
	# points(all.res[, "partitions"], all.res[, "aic"], col = col[["aic"]], type = "l");
	# points(all.res[, "partitions"], all.res[, "aic"], col = col[["aic"]], pch = 21, bg = "white");

	# legend("topleft", c("aicc", "aic"), pch = 21, pt.bg = "white", lty = 1, col = c(col[["aicc"]], col[["aic"]]), inset = 0.05, cex = 0.75, bty = "n");
# }




# ## Consider removing previously-fit rate shifts
# #.back.step.medusa <- function (currentModel, z, step, model, fixPar, criterion) {
.back.step.medusa <- function (currentModel, z, step, model, criterion) {
## As a first step, only consider removing entire shifts. Later deal with individual parameters.
	z.opt <- z;
	bestModel <- currentModel;
	bestScore <- as.numeric(bestModel[criterion]);
	allDeletedShifts <- NULL;
	bestRemoved <- NULL;
	improve <- T;
	
	while (improve) { # may be possible to remove > 1 previously fit shift
		allDeletedShifts <- c(allDeletedShifts, bestRemoved);
		currentModel <- bestModel;
		z <- z.opt;
		cuts <- bestModel$cut.at;
		nodes <- bestModel$split.at;
		pars <- bestModel$par;
		numModels <- length(bestModel$par)/2;
		improve <- F;
		
		if (numModels > 2) {
			for (i in 2:(numModels - 1)) { # don't waste time removing last shift
				fitModel <- currentModel;
				obj <- .dissolve.split.medusa(z, cut = cuts[i], node = nodes[i], aff = i);
				aff <- obj$affected;
				z.temp <- obj$z[obj$z[,"partition"] == aff,,drop=FALSE];
				
		## set par to weighted mean of 2 affected partitions
			## updated to reflect total path length rather than number of tips
			## really only influences weighted parameter starting values
			
				weights <- c(sum(z[which(z[,"partition"] == aff),"t.len"]), sum(z[which(z[,"partition"] == i),"t.len"]));
				sp <- c(weighted.mean(pars[c(aff,i),1], weights), weighted.mean(pars[c(aff,i),2], weights, na.rm=T));
				#fit <- getOptimalModelFlavour(z=z.temp, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
				fit <- .get.optimal.model.flavour(z=z.temp, sp=sp, model=model, criterion=criterion);
				#.fit.partition.medusa
				
				#.fit.partition.medusa(partition = 2, z = z.bd.stem, sp = c(initialR, initialE), model = "bd");
				
		## Update fit values
				fitModel$par[aff,] <- fit$par;
				fitModel$par <- fitModel$par[-i,];
				fitModel$lnLik.part[aff] <- fit$lnLik;
				fitModel$lnLik.part <- fitModel$lnLik.part[-i];
				fitModel$lnLik <- sum(fitModel$lnLik.part);
				model.fit <- .calculate.modelfit.medusa(fit=fitModel, z=z);
				fitModel$aic <- model.fit$aic;
				fitModel$aicc <- model.fit$aicc;
				fitModel$num.par <- model.fit$k;
				
				if (fitModel[criterion] < bestScore) {
					#cat("Found a better fit by removing a model (split at node #", fitModel$split.at[i], ")!\n", sep="");
					fitModel$split.at <- fitModel$split.at[-i];
					fitModel$model[aff] <- fit$model;
					fitModel$model <- fitModel$model[-i];
					fitModel$cut.at <- fitModel$cut.at[-i];
					bestModel <- fitModel;
					bestScore <- as.numeric(fitModel[criterion]);
					z.opt <- .updateZ(z=obj$z, deletedPart=i);
					bestRemoved <- nodes[i];
					improve <- T;
				}
			}
			if (improve) {step <- rbind(step, c("remove", bestRemoved));}
		}
	}
	return(list(fit=bestModel, z=z.opt, step=step, remove=bestRemoved));
}


## Remove previously-fit rate shift
.dissolve.split.medusa <- function (z, cut, node, aff) {
## Grab ancestral branch partition membership
	anc <- z[which(z[,"dec"] == node)];
	root <- min(z[,"anc"]);
	tag <- NULL;
	
	if (cut == "node") {
		tag <- as.numeric(z[which(z[,"dec"] == node),"partition"]);
	} else if (cut == "stem" && anc > root) {
		tag <- as.numeric(z[which(z[,"dec"] == anc),"partition"]);
	} else if (cut == "stem" && anc == root) { # need to take other side of root
		dec <- z[which(z[,"anc"] == root),"dec"];
		tag <- as.numeric(z[which(z[,"dec"] == dec[which(dec != node)]),"partition"]); # ug. li.
	}
	
	idx <- which(z[,"partition"] == aff);
	z[idx,"partition"] <- tag;
	
	return(list(z=z, affected=tag));
}

## Only print if model improves AIC score
.print.removed.shifts <- function (remove) {
	for (i in 1:length(remove)) {
		cat("  Removing shift at node #", remove[i], "\n", sep="");
	}
}

## reduce partition IDs to reflect dissolved split
.updateZ <- function (z, deletedPart) {
	idx <- z[,"partition"] > deletedPart;
	z[idx,"partition"] <- z[idx,"partition"] - 1;
	return(z);
}

## passed-in parameters parm are MLEs stored in a matrix
.get.profile.likelihoods <- function (res, crit=1.92) {
	parm <- res$model$par;
	models <- res$model$model;
	z <- res$model$z;
	
	prof.par <- matrix(nrow=length(parm[,1]), ncol=4);
	colnames(prof.par) <- c("r.low", "r.high", "eps.low", "eps.high")
	inc <- 0.05;
	
	# get this out of for-loop format
	for (i in 1:length(parm[,1])) {
		#cat("Model",i,"\n")
		model <- models[i]
		sp <- as.numeric(parm[i,]);
		new.part <- z[z[,"partition"] == i,,drop=FALSE];
		lik <- .lik.partition.medusa(partition=new.part, model=model);
		
		if (model == "yule") {
			par <- sp[1];
			maxLik <- lik(par); if (maxLik == -Inf) maxLik <- 0; # correct for -Inf at boundary lambda == 0
			
			thprof.parhold <- function (x) lik(x) - maxLik + crit; # find roots on either side of maxLik
			
	## need intelligent bounds
			if (par != 0) {
				low.bound <- par - par/2;
				up.bound <- par + par/2;
			} else {
				low.bound <- par;
				up.bound <- par + inc/2;
			}
			if (low.bound != 0) {
				while (thprof.parhold(low.bound) > 0) {
					low.bound <- low.bound - inc;
				}
			}
			while (thprof.parhold(up.bound) > 0) {
				up.bound <- up.bound + inc;
			}
			
			if (low.bound <= 0) {
				prof.par[i,1] <- 0;
			} else {
				prof.par[i,1] <- uniroot(thprof.parhold, lower=low.bound, upper=par)$root;
			}
			if (par == 0) par <- 1e-10; # avoid -Inf at boundary
			prof.par[i,2] <- uniroot(thprof.parhold, lower=par, upper=up.bound)$root;
			
		} else if (model == "bd") {
			par1 <- sp[1]; par2 <- sp[2];
			maxLik <- lik(sp);
			
	## first, r
			thprof.parholdR <- function (x) lik(c(x, par2)) - maxLik + crit;
			
			low.bound <- par1 - par1/2;
			up.bound <- par1 + par1/2;
			
			while (thprof.parholdR(low.bound) > 0) {
				low.bound <- low.bound - inc;
			}
			while (thprof.parholdR(up.bound) > 0) {
				up.bound <- up.bound + inc;
			}
			if (low.bound <= 0) low.bound <- 0;
			
			prof.par[i,1] <- uniroot(thprof.parholdR, lower=low.bound, upper=par1)$root;
			prof.par[i,2] <- uniroot(thprof.parholdR, lower=par1, upper=up.bound)$root;
			
	## now, epsilon
			thprof.parholdE <- function (x) lik(c(par1, x)) - maxLik + crit;
			
			low.bound <- par2 - par2/2;
			up.bound <- par2 + par2/2;
			
			while (thprof.parholdE(low.bound) > 0) {
				low.bound <- low.bound - inc;
			}
			while (thprof.parholdE(up.bound) > 0) {
				up.bound <- up.bound + inc;
			}
			if (low.bound < 0) low.bound <- 0;
			if (up.bound > 1) up.bound <- 1;
			
			if (low.bound == 0) {
				prof.par[i,3] <- 0;
			} else {
				prof.par[i,3] <- uniroot(thprof.parholdE, lower=0, upper=par2)$root;
			}
			if (up.bound == 1) {
				prof.par[i,4] <- 1;
			} else {
				prof.par[i,4] <- uniroot(thprof.parholdE, lower=par2, upper=up.bound)$root;
			}
			if (prof.par[i,3] == par2) prof.par[i,3] <- 0; # precision problem?!? check optimization ***
		}
	}
	# remove columns that are not used
	idx <- which(apply(prof.par, MARGIN=2, function (x) all(is.na(x))));
	if (length(idx) > 0) {
		prof.par <- prof.par[,-idx];
	}
	prof.par <- as.data.frame(round(prof.par, digits=7));
	return(prof.par);
}


## Create a plot of model-fit vs. model-size
#plotModelFit
#.plot.modelfit.medusa <- function (all.res, ...)
#{
#       ylim <- c(min(all.res[,"aic"],all.res[,"aicc"]), max(all.res[,"aic"],all.res[,"aicc"]));
#       plot(all.res[,"partitions"],all.res[,"aicc"], xlab="number of piecewise models", ylab="model fit", ylim=ylim, type="l", col="blue", xaxt="n", ...);
#   axis(1, at=all.res[,"partitions"], labels=all.res[,"partitions"], tcl=NA)
#       points(all.res[,"partitions"],all.res[,"aicc"], col="blue", pch=21, bg="white");
#       points(all.res[,"partitions"],all.res[,"aic"], col="black", type="l");
#       points(all.res[,"partitions"],all.res[,"aic"], col="black", pch=21, bg="white");

#       legend("topleft", c("aicc","aic"), pch=21, pt.bg="white", lty=1, col=c("blue", "black"), inset = .05, cex=0.75, bty="n"); # 'bottomright' also works
#}

## Takes in a summary from summarizeTurboMEDUSA
## treeParameters <- list(mm=mm, break.pts=break.pts, phy=phy, z=z)
#plotPrettyTree <- function (treeParameters, time=TRUE, node.labels=FALSE, margin=FALSE, cex=0.5, label.offset=0, font=3, color.tip.label=FALSE, ...)
#{
#       mm <- treeParameters$mm;
#       break.pts <- treeParameters$break.pts;
#       phy <- treeParameters$phy;
#       z <- treeParameters$z;
#       colour <- NULL;

#       dev.new();
#
#       if (color.tip.label)
#       {
#               for (i in 1:length(phy$tip.label))
#               {
#                       colour[i] <- as.integer(z[which(z[,"dec"] == i),"partition"]);
#               }
#       }
#       if (time) {margin=TRUE;}
#       plot.phylo(phy, edge.color=z[mm,"partition"], no.margin=!margin, cex=cex, label.offset=label.offset, tip.color=colour, ...);
#       if (time)
#       {
#               axisPhylo(cex.axis=0.75);
#               mtext("Divergence Time (MYA)", at=(max(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx)*0.5), side = 1, line = 2, cex=0.75);
#       }
#       if (node.labels)
#       {
#               for (i in  1:length(break.pts))
#               {
#                       nodelabels(i, node= break.pts[i], frame = "c", font = 1, cex=0.5);
#               }
#       }
#}


## Get b and d values from r (b-d) and epsilson (d/b)
## Used in previous version of program; now in terms of r and epsilon
## Possibly of use to users wishing to translate results
#get.b.d <- function (r, epsilon)
#{
#       b <- r/(1-epsilon);
#       d <- b-r;   # Alternatively: d <- eps*r/(1-eps)
#       return(list(b=b, d=d));
#}

## Print out tree with ape-style node-numbering
## Possibly of interest for users to identify numbers of node(s) off interest
## If this is the case, make sure to pass in pruned tree
#plotNN <- function (phy, time=TRUE, margin=TRUE, label.offset=0.5, cex=0.5, ...)
#{
#       phy$node.label <- (length(phy$tip.label) + 1):max(phy$edge);
#       plot.phylo(phy, show.node.label=TRUE, no.margin=!margin, label.offset=label.offset, cex=cex, ...);
#       if (time && !margin) {cat("Cannot plot time axis without a margin.\n");}
#       else if (time && margin) {axisPhylo(cex.axis=0.75)};
#}