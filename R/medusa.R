#takes a phylogeny where the tips of the phylogeny may represent stems of unresolved clades and a list of taxonomic richness for each tip clade and fits a series of birth-death models to each branch in the tree.   Default is estimate d/b. Model limit is the upper limit on splits to be tried. If AIC scores are still improving substanitally after 20 splits, this limit should be increased. 
runMedusa <- function (phy, richness, estimateExtinction=T, modelLimit=20, cutAtStem=T, startR=0.05, startE=0.5, ...) 
{

	if(dim(richness)[2]==2) {
		nr<-dim(richness)[1]
		richness<-data.frame(richness[,1], rep(0, nr), richness[,2])
	}
    phy$node.label <- NULL
    
    #reset model limit if tree has fewer than 20 internal branches
	N<-length(phy$tip.label)
	if(modelLimit>(2*N-2)) modelLimit=2*N-2
   
  #holds parameter estimates  
    allRes<-matrix(nrow=modelLimit+1, ncol=5)
    
    root <- max(phy$edge) - phy$Nnode + 1
    node.list <- 1:max(phy$edge)
    node.list <- node.list[node.list != root]

#First pass estimates a birth-death model that integrates phylogenetic and taxonomic likelihoods described in Rabosky et al., 2007      	  
	  
	res <- list()

	baseModel<-fitDiversification(phy, richness, estimateExtinction=T)

      allRes[1,]<-c(0, baseModel$LH, baseModel$np, baseModel$aic, baseModel$aicc)
      
	#Now split the tree on all branches. 

      for (i in 1:length(node.list)) {
        r1 <- NULL
        r2 <- NULL
        z1 <- NULL
        z2 <- NULL
        z <- NULL
        z <- splitEdgeMatrixGeiger(phy, node.list[i], richness, cutAtStem)
        z1 <- z[z[, 7] == 1, ]
        z2 <- z[z[, 7] == 2, ]
        
        if(length(dim(z1))==0) z1<-rbind(z1)
        if(length(dim(z2))==0) z2<-rbind(z2)

	#Fit BD model to the two subclades of the tree        

        r1 <- getDiversificationModel(z1,  estimateExtinction, startR, startE)
        r2 <- getDiversificationModel(z2, estimateExtinction, startR, startE)

   	 	if(estimateExtinction) np<-5 else np<-3   
    
   
        k<-2*nrow(richness)-1

        res$node[i] <- node.list[i]
        res$LH[i] <- r1$LH + r2$LH
        res$aic[i] <- (-2 * res$LH[i]) + 2*np
        res$aicc[i]<-res$aic[i]+2*np*(np+1)/(k-np-1)

        if(estimateExtinction) {
        	eps1<-r1$par[2]
        	eps2<-r2$par[2] 
 
        } else {
        	eps1=0
            eps2=0
        }
        
        
        res$r1[i] <- r1$par[1]
        res$e1[i] <- eps1
        res$LH1[i] <- r1$LH
        
        res$r2[i] <- r2$par[1]
        res$e2[i] <- eps2
        res$LH2[i] <- r2$LH
        
        res$np[i]<-np
        

      }

      bestModel<-which(res$LH==max(res$LH))[1]
      allRes[2,]<-c(res$node[bestModel], res$LH[bestModel], res$np[bestModel], res$aic[bestModel], res$aicc[bestModel])
      z<-splitEdgeMatrixGeiger(phy, node.list[bestModel], richness, cutAtStem)
      
    
    for(j in 2:modelLimit) {
    	
      res <- list()
      for (i in 1:length(node.list)) {
 
        zNew <- resplitEdgeMatrixGeiger(z, phy, node.list[i], cutAtStem)# resplit edge matrix adds a new split to an already-split tree
        LH=0
        np=0
        for(k in 1:max(zNew[,7])){
        	zPart <- zNew[zNew[, 7] == k, ]
        	if(length(dim(zPart))==0) zPart<-rbind(zPart)
        	if(nrow(zPart)!=0) {
        		rPart <- getDiversificationModel(zPart, estimateExtinction, startR, startE)

        		LH<-LH+rPart$LH
        		if(estimateExtinction) np<-np+3 else np<-np+2
        	}
        	
		}
		np<-np-1 # otherwise you charge one extra for the background
    
   
             k<-2*nrow(richness)-1


        res$node[i] <- node.list[i]
        res$LH[i] <- LH
        res$aic[i] <- (-2 * res$LH[i]) + 2*np
        res$aicc[i]<-res$aic[i]+2*np*(np+1)/(k-np-1)

        res$np[i] <- np

      }

      bestModel<-which(res$LH==max(res$LH))[1]
      allRes[j+1,]<-c(res$node[bestModel], res$LH[bestModel], res$np[bestModel], res$aic[bestModel], res$aicc[bestModel])
      
      z <- resplitEdgeMatrixGeiger(z, phy, node.list[bestModel], cutAtStem)
      cat(j, "\n")
      
    }
    
    return(allRes)
}

# Adds a new split to an already-split edge matrix
resplitEdgeMatrixGeiger <- function (z, phy, node, cutAtStem=T) 
{
	newTag<-max(z[,7])+1
    rootnode <- length(phy$tip.label) + 1

    if (node >= rootnode) {
        node.desc <- node
        pos <- 1
        
        if(cutAtStem) {
        	row1<-which(z[, 2] == node.desc[1])
        	row2<-which(z[, 1] == node.desc[1])
 			row<-c(row1, row2)
 		} else {
 			row<-which(z[, 2] == node.desc[1])
 			}       
        base<-min(z[row,7])
        ok<- z[row,7]==base 
        if(sum(ok)>0)
               z[row[ok],7] <- newTag
           

        while (pos != (length(node.desc) + 1)) {
            temp <- node.sons(phy, node.desc[pos])
            temp <- temp[temp > rootnode]
            if(length(temp)!=0) {
             for (k in 1:length(temp)) {
             	row<-which(z[, 1] == temp[k])
                ok<- z[row,7]==base 
                if(sum(ok)>0)
                	z[row[ok],7] <- newTag
            	}
             node.desc <- c(node.desc, temp)
            }
            pos <- pos + 1
        }
    }
    else if (node > 0) 
        z[z[, 2] == node,7] <- newTag

    return(z)
}


#This returns b and d parameters from composites r (net diversification) and epsilson (death over birth rate)
getBD<-function(r, eps) {
	d<-eps*r/(1-eps)
	b<-r+d
	return(list(b=b, d=d))
	}
	
#This takes a phylogeny and a list of branches where rates shift and returns a list that describes the shift points and new rates.

getFullSplitModel<-function(phy, estimateExtinction=T, breakList, richness, cutAtStem=T){

	res<-list()
	phy$node.label <- NULL
    root <- max(phy$edge) - phy$Nnode + 1
    node.list <- 1:max(phy$edge)
    node.list <- node.list[node.list != root]
    x <- branching.times(phy)

	# First break

    z <- splitEdgeMatrixGeiger(phy, breakList[1], richness, cutAtStem=cutAtStem)
	if(length(breakList)>1)
	  for(i in 2:length(breakList)) {
		z <- resplitEdgeMatrixGeiger(z, phy, breakList[i], cutAtStem=cutAtStem)

		}
	
	 for(k in 1:max(z[,7])){
        	zPart <- z[z[, 7] == k, ]
        	if(length(dim(zPart))==0) zPart<-rbind(zPart)

        	if(nrow(zPart)!=0) {
        		rPart <- getDiversificationModel(zPart, estimateExtinction)
        		res[[k]]<-rPart
        		int<-is.na(zPart[,6])
        		if(sum(int)==0)  {
        			nn<-phy$tip.label[zPart[,2]]
        			
        	        res[[k]]$taxa<-nn

        		} else {
        			tt<-zPart[!int, 2]
        			res[[k]]$taxa<-phy$tip.label[tt]

        			}
        
        	}
        	
		}
		list(res=res, z=z)
}	


fitDiversification <- function (phy, richness, estimateExtinction=T, ...) 
{
    z <- splitEdgeMatrixGeiger(phy, NA, richness, cutAtStem=T)
    r1 <- getDiversificationModel(z, estimateExtinction, ...)
    if(estimateExtinction) np<-2 else np<-1    
    res <- list()
    res$LH <- r1$LH
    res$aic <- (-2 * r1$LH) + 2*np
    n<-nrow(z)
    res$np<-np
    res$k<-n
    res$aicc<-res$aic+2*np*(np+1)/(n-np-1)
    if(estimateExtinction) eps<-r1$par[2] else eps=0
    
    res$r <- r1$par[1]
    res$eps <- eps
    return(res)
}

getDiversificationModel <- function (z, estimateExtinction=T,     startR=0.05, startE=0.5) 
{
    isInt <- is.na(z[,6])
    isTerm <- !isInt
    int<-z[isInt,]
    term<-z[isTerm,]
    
    if(length(dim(int))==0) int<-rbind(int)
    if(length(dim(term))==0) term<-rbind(term)

    nint <- sum(isInt)
    nterm <- sum(isTerm)
    
    lalphaF<-function(b, d, aa, tt) {
    	num<-log(d)+log(exp((b-d)*tt)-1)
    	den<-log(b*(exp((b-d)*tt))-d)     	
    	res<-aa*(num-den)
    	res
    	}	
    
    lbetaF <- function(b, d, aa, tt) {
    	res<-lalphaF(b, d, aa, tt)+log(b)-log(d)
    	res
    }	
    
    lbetaPB <- function(b, aa, tt) {
    	res<-log((exp(b*tt)-1)/(exp(b*tt)))
    }


    Lfunc_comb <- function(r, eps) {
    	if(r<0 | eps<=0 | eps>=1) return(-Inf)

    	bd<-getBD(r, eps)
    	b<-bd$b
    	d<-bd$d
    	
    	
    	if(nint==0) {lint=0} else {
    		branchLengths<-int[,3]-int[,4]
    		lint<-nint * 
            log(r) - r * sum(branchLengths) - sum(log(1 - (eps * 
            exp(-r * int[, 3]))))
        }
    	if(nterm==0) {lterm=0} else {
    		a<-term[,"startR"]
    		n<-term[,"endR"]
    		timeInterval<-term[,"startTime"]-term[,"endTime"]
    		endj<-pmin(a, n)
    		sum<-0
    		lnl<-numeric(nterm)
    		for(i in 1:nterm) {
    			lnxx<-numeric(length=endj[i])
    			for(j in 1:endj[i]) { #using likelihoods from Foote et al. 1999, a correction of Raup
    				logAlpha<-lalphaF(b, d, a[i], timeInterval[i])
    				logBeta<-lbetaF(b, d, a[i], timeInterval[i])
    				if(logBeta>0) logBeta=0
    				s1<-lchoose(a[i],j)+lchoose(n[i]-1,j-1)
    				
    				if(logAlpha==-Inf) s2<-0 else s2<-(a[i]-j)*logAlpha
    				s3<-log(((1-exp(logAlpha))*(1-exp(logBeta)))^j)
    				s4<-(n[i]-1)*logBeta
    				s5<-log(1-exp(logAlpha)) # Conditioning on survival to the present
       				lnxx[j]<-s1+s2+s3+s4-s5

    			}
    			lnl[i]<-logspace_sum(lnxx)
    		}
    		lterm<-sum(lnl)
    	}
    	#return(lterm)

    	return(lint+lterm)
	}
	
	Lfunc_comb_pb <- function(r) {
    	if(r<0) return(-Inf)

    	b<-r
    	d<-0
    	eps<-0
    	
    	if(nint==0) {lint=0} else {
    		branchLengths<-int[,3]-int[,4]
    		lint<-nint * 
            log(r) - r * sum(branchLengths) - sum(log(1 - (eps * 
            exp(-r * int[, 3]))))
        }
    	if(nterm==0) {lterm=0} else {
    		a<-term[,"startR"]
    		n<-term[,"endR"]
    		timeInterval<-term[,"startTime"]-term[,"endTime"]
    		endj<-pmin(a, n)
    		sum<-0
    		lnl<-numeric(nterm)
    		for(i in 1:nterm) {
    			lnxx<-numeric(length=endj[i])
    			for(j in 1:endj[i]) { #using likelihoods from Foote et al. 1999, a correction of Raup
    				logAlpha<- -Inf
    				logBeta<-lbetaPB(b, a[i], timeInterval[i])
    				s1<-lchoose(a[i],j)+lchoose(n[i]-1,j-1)
    				
    				if(logAlpha==-Inf) s2<-0 else s2<-(a[i]-j)*logAlpha
    				s3<-log(((1-exp(logAlpha))*(1-exp(logBeta)))^j)
    				s4<-(n[i]-1)*logBeta
    				s5<-log(1-exp(logAlpha)) # Conditioning on survival to the present
       				lnxx[j]<-s1+s2+s3+s4-s5

    			}
    			lnl[i]<-logspace_sum(lnxx)
    		}
    		lterm<-sum(lnl)
    	}
    	
    	return(lint+lterm)
	}

	
    res <- list()
    

    
    if (estimateExtinction == TRUE) {
     	foo <- function(x) 
     		- Lfunc_comb(r=exp(x[1]), eps=exp(x[2]))
		sp<-log(c(startR, startE))
		o<-optim(foo, par=sp, method="N")
		res$LH <- -o$value
    	res$par <- exp(o$par)

    } else {
    	foo <- function(x) 
     		- Lfunc_comb_pb(r=exp(x[1]))
		sp<-log(startR)
		o<-optimize(foo, interval=c(-100, 1))
		res$LH <- -o$objective
    	res$par <- exp(o$minimum)
    	}

	    



    return(res)
}


plotDiversificationSurface <- function (phy, nPoints=10, rInterval=c(0.00001, 0.3), eInterval=NULL, logTransform=T) 
{
    z <- splitEdgeMatrixGeiger(phy, phy$Nnode)
	rootnode = (length(phy$tip.label) + 1)
 	int <- z[z[, 2] > rootnode, ]
    term <- z[z[, 2] < rootnode, ]
    nint <- nrow(int)
    nterm <- nrow(term)
    betaF <- function(r, eps, t1) {
    	if(r<0 | eps<0 | eps>=1) return(-Inf)
        xf <- (exp(r * t1) - 1)/(exp(r * t1) - eps)
        xf
    }

    Lfunc_comb <- function(r, eps) {
    	if(r<0 | eps<0 | eps>=1) return(-Inf)

        (sum(log(1 - betaF(r, eps, term[1:nterm, 4]))) + sum((term[1:nterm, 
            5] - 1) * log(betaF(r, eps, term[1:nterm, 4]))) + nint * 
            log(r) - r * sum(int[1:nint, 4]) - sum(log(1 - (eps * 
            exp(-r * int[1:nint, 3])))))
    }
	
	if(logTransform) {
			x<-exp(seq(from=log(rInterval[1]), to=log(rInterval[2]), length.out=nPoints))
	} else x<-seq(from=rInterval[1], to=rInterval[2], length.out=nPoints)
	
	if(is.null(eInterval)) {
		lnl<-numeric(nPoints)
		for(i in 1:nPoints)
			lnl[i]<-Lfunc_comb(r=x[i], eps=0)
		plot(x, lnl, type="l")	
	} else {
		if(logTransform) {
			y<-exp(seq(from=log(eInterval[1]), to=log(eInterval[2]), length.out=nPoints))
		} else y<-seq(from=eInterval[1], to=eInterval[2], length.out=nPoints)
		lnl<-matrix(nrow=nPoints, ncol=nPoints)
		for(i in 1:nPoints)
			for(j in 1:nPoints)
				lnl[i,j]<-Lfunc_comb(r=x[i], eps=y[j])
		mm<-max(lnl)		
		contour(lnl, levels=c(mm-1, mm-2, mm-3, mm-4, mm-5, mm-10, mm-100), labels=c(1, 2, 3, 4, 5, 10, 100), axes=F, xlab="b-d", ylab="d/b")
		tics<-floor(c(1, nPoints/4, nPoints/2, nPoints*3/4, nPoints))
		axis(1, at=c(0, 0.25, 0.5, 0.75, 1), labels=round(x[tics], 3))
		axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=round(y[tics], 3))
	}
}




fitSplitModel <- function (phy, estimateExtinction=T,     startR=0.05, startE=0.5) 
{
    phy$node.label <- NULL
    root <- max(phy$edge) - phy$Nnode + 1
    node.list <- 1:max(phy$edge)
    node.list <- node.list[node.list != root]
    x <- branching.times(phy)
    res <- list()
    for (i in 1:length(node.list)) {
        r1 <- NULL
        r2 <- NULL
        z1 <- NULL
        z2 <- NULL
        z <- NULL
        z <- splitEdgeMatrixGeiger(phy, node.list[i])
        z1 <- z[z[, 6] == 1, ]
        z2 <- z[z[, 6] == 2, ]
        
        r1 <- getDiversificationModel(z1, estimateExtinction, startR, startE)
        r2 <- getDiversificationModel(z2, estimateExtinction, startR, startE)

   	 	if(estimateExtinction) np<-4 else np<-2    
    
   
        

        res$node[i] <- node.list[i]
        res$LH[i] <- r1$LH + r2$LH
        res$aic[i] <- (-2 * res$LH[i]) + 2*np*2
        
        if(estimateExtinction) {
        	eps1<-r1$par[2]
        	eps2<-r2$par[2] 
 
        } else {
        	eps1=0
            eps2=0
        }
        
        
        res$r1[i] <- r1$par[1]
        res$e1[i] <- eps1
        res$LH1[i] <- r1$LH
        
        res$r2[i] <- r2$par[1]
        res$e2[i] <- eps2
        res$LH2[i] <- r2$LH

    }
    res <- as.data.frame(res)
    return(res)
}




splitEdgeMatrixGeiger <- function (phy, node, richness, cutAtStem=T) 
{
    bt <- branching.times(phy)
    rootnode <- length(phy$tip.label) + 1
    
    #m1<-match(richness[,1], phy$tip.label)
    #m2<-match(phy$edge[,1], names(bt))

    interior<-phy$edge[,2] %in% phy$edge[,1]
 	tips<-!interior 
 	
 	z<-matrix(ncol=7, nrow=sum(interior)+nrow(richness))
    colnames(z)<-c("anc", "dec", "startTime", "endTime", "startR", "endR", "Group")
  
 	for(i in 1:sum(interior)) {
 		anc<-phy$edge[interior,1][i]
 		dec<-phy$edge[interior,2][i]
 		startTime<-bt[names(bt)==anc]
 		endTime<-startTime-phy$edge.length[interior][i]
 		sr<-1
 		er<-NA
 		group<-NA
 		z[i,]<-c(anc, dec, startTime, endTime, sr, er, group)
 		} 
  
  	m1<-match(richness[,1], phy$tip.label)

  
	for(i in 1:nrow(richness)) {
		wt<-phy$tip.label[m1][i]
		we<-which(phy$edge[,2]==m1[i])
		anc<-phy$edge[we,1]
 		dec<-phy$edge[we,2]
 		startTime<-bt[names(bt)==anc]
 		endTime<-richness[i,2]
 		sr <- 1
 		er<-richness[i,3]
 		group<-NA
 		z[i+sum(interior),]<-c(anc, dec, startTime, endTime, sr, er, group)
		}
	 
	for(i in 1:max(z[,2])) {
		nn<-which(z[,2]==i)
		o<-order(z[nn,4], decreasing=T)
		if(length(nn)>1) for(j in 2:length(nn)) {
			z[nn,][o[j],3]<-z[nn,][o[j-1],4]
			z[nn,][o[j],5]<-z[nn,][o[j-1],6]
			}
		}
	if(is.na(node)) {
		z[,7]<-1
	} else {
		setZ<-function(z, n) {
			if(n %in% z[,2]) {
       			ok<-z[,2]==n
       			z[ok,7]<-2
				ns<-node.sons(phy, n)
				if(length(ns)!=0) for(k in 1:length(ns))
					z=setZ(z, ns[k])
			} else if(n %in% z[,1]) z[,7]<-2
			z
		}
		
		z<-setZ(z, node)
        
    	z[is.na(z[,7]),7]<-1
    }            
    
    if(!cutAtStem) {
    	row<-which(z[, 1] == node)
    	if(length(row)>0)
			z[row,7]<-1
    	}
       
    return(z)
}
##########################



summaryMedusa<-function(phy, richness, out, cutoff=4, plotTree=T, useCorrection=F, cutAtStem=T) {
	breaks<-numeric()
	i=2
	if(useCorrection) col=5 else col=4
	while(1) {
		if((out[i-1,col]-out[i,col])<cutoff) break;
		breaks[i-1]<-out[i,1]
		i<-i+1
		}
	rr<-getFullSplitModel(phy, breakList=breaks, richness=richness, cutAtStem=cutAtStem)
	
	if(plotTree) {
		mm<-match(phy$edge[,2], rr$z[,2])
		plot(phy, edge.color=rr$z[mm,7])
		}
	
	
	rr$res	
	}
