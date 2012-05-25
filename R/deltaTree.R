############################################################
###############TREE TRANSFORMATIONS#########################
############################################################
## a few functions modified 11 dec 07 to accommodate ultrametric trees; see comments below

###LAMBDA###
## modified 11 dec 07
lambdaTree <- function(phy, lambda)
{
	originalDiag = diag(vcv.phylo(phy)) ## added to fix the final edge length problem in non-ultrametric tree
	original.rtt <- max(diag(vcv.phylo(phy))) ## generalizes to a non-ultrametric tree
	  ## following line assumes ultrametricity
	  ## original.rtt <- max(branching.times(phy))
	ltree <- phy
	ltree.old <- new2old.phylo(phy) #needed due to change in ape tree format
	ltree$edge.length <- ltree$edge.length * lambda #shortens internal branches in proportion to lambda
	  ## The following assignment of terminal branch fails when the tree is non-ultrametric
	  ## t <- original.rtt * (1-lambda) #needed to rescale tree back to its original root-to-tip length
	which(ltree.old$edge[,2] > 0) -> terminal.edges #identifies termal edges to be extended for recovery of original rtt length; based on old tree format in which nodes are coded as negative values
	terminalEdgesOrdered <- match(match(labels(originalDiag), ltree$tip.label, nomatch = F), ltree$edge[,2])
	ltree$edge.length[terminalEdgesOrdered] <- 
	ltree$edge.length[terminalEdgesOrdered] + (originalDiag - (originalDiag * lambda)) ## generalizes to a non-ultrametric tree
	  ## The following didn't work on a non-ultrametric tree
	  ## ltree$edge.length[terminal.edges] + t -> ltree$edge.length[terminal.edges] -- extends terminal edges
	estMinLambda <- 
	  ((min(phy$edge.length[terminalEdgesOrdered]) 
	  + originalDiag[phy$edge[,2][match(min(phy$edge.length[terminalEdgesOrdered]), phy$edge.length)]]) 
	  / originalDiag[phy$edge[,2][match(min(phy$edge.length[terminalEdgesOrdered]), phy$edge.length)]] )
	if(any(ltree$edge.length < 0)) message("Lambda values too large imply negative terminal branch lengths. \nLambda should not exceed an estimated ", round(estMinLambda, 3), " on this tree.")
	return(ltree)
}

###DELTA###
## checked 11 dec 07; works fine
## Modified Jan 2. 08 by Steve Kembel to work when trees have internal branch labels

deltaTree <- function (phy, delta, rescale = T) {    tmp <- as.numeric(phy$edge)    times <- branching.times(phy)    original.rtt.length <- max(times)    times = max(times) - times    max(times)    tips <- length(phy$tip.label)    res <- phy    for (i in 1:length(phy$edge.length)) {        bl <- phy$edge.length[i]        age = times[phy$edge[i, 1] - tips]        res$edge.length[i] <- (age + bl)^delta - age^delta    }    if (rescale == T)         res <- rescaleTree(res, original.rtt.length)    res}

###TWORATE#####

tworateTree<-function(phy, breakPoint, endRate) 
{
	times<-branching.times(phy)	
	for(i in 1:length(phy$edge.length)) {
		bl<-phy$edge.length[i]
		age=times[which(names(times)==phy$edge[i,1])] #gets tip to node length
		if((age-bl)<breakPoint) #identifies branches that are on the tip side of the break-point (i.e., young)
			phy$edge.length[i]<-(age-min(age, breakPoint))*1+(min(age, breakPoint)-(age-bl))*endRate #If an edge is entirely to the right of the breakpoint it is simply multipled by f.  However, if the edge extends across the breakpoint, we want to leave the part to the left of the BP unchanged (this is the first part of this equation, which multiplies this part by 1).  We then want to multiply the part of the edge that is to the right of the BP by f (this is the second part of this line).
		}
	phy
}

###LINEAR CHANGE###
linearchangeTree<-function(phy, endRate=NULL, slope=NULL)
{
	
    if(is.null(slope)&&is.null(endRate))
    	stop("Must supply either endRate or slope")
    	
	times<-branching.times(phy)
	names(times)<-(as.numeric(names(times)))
	rootdepth=max(times)
	
	
	
	if(is.null(slope)) {
		slope=(endRate-1)/rootdepth
		}
	

	getrate<-function(begin, end)
	{
		br<-1+begin*slope
		er<-1+end*slope
		if(br>0 & er>0) return((br+er)/2)
		if(br<0 & er<0) return(0)
		ii<- -1/slope
		return(br*(ii-begin)/(2*(end-begin)))
	}

	for(i in 1:length(phy$edge.length)) {
		bl<-phy$edge.length[i]
		age=rootdepth-times[which(names(times)==phy$edge[i,1])]
		end=age+bl
	
		phy$edge.length[i]<-phy$edge.length[i]*getrate(age, end)

		

	}
	phy	
}

###RESCALING###
## modified 11 dec 07
rescaleTree<-function(phy, totalDepth)
{
	## following line added to accommodate ultrametric trees
	d <- max(diag(vcv.phylo(phy)))
	  ## assumes ultrametricity
	  ## d<-max(branching.times(phy))
	phy$edge.length<-(phy$edge.length/d)*totalDepth
	phy
}

exponentialchangeTree<-function (phy, endRate=NULL, a=NULL) 
{
	

    if(is.null(a)&&is.null(endRate))
    	stop("Must supply either endRate or a")
    	
   	times <- branching.times(phy)
    d<-max(times)
    if(is.null(a))
 	   a<-log(endRate)/d
 	if(a==0) return(phy)   
    names(times) <- (as.numeric(names(times)))
    for (i in 1:length(phy$edge.length)) {
        bl <- phy$edge.length[i]
        age = times[which(names(times) == phy$edge[i, 1])]
        t1 = max(times) - age
        t2 = t1+bl
        phy$edge.length[i] = (exp(a*t2)-exp(a*t1))/(a)
    }
    phy
}

speciationalTree<-function(phy)
{
	phy$edge.length[]=1
	phy	
}

ouTree<-function(phy, alpha) {
	times <- branching.times(phy)
    names(times) <- (as.numeric(names(times)))
    Tmax<-times[1]
    phy2<-phy
    for (i in 1:length(phy$edge.length)) {
        bl <- phy$edge.length[i]
        age = times[which(names(times) == phy$edge[i, 1])]
        t1 = max(times) - age
        t2 = t1+bl
        phy2$edge.length[i] = (1/(2*alpha))*exp(-2*alpha * (Tmax-t2)) * (1 - exp(-2 * alpha * t2)) - 
        						(1/(2*alpha))*exp(-2*alpha * (Tmax-t1)) * (1 - exp(-2 * alpha * t1))
    }
    phy2
	
}

kappaTree<-function(phy, kappa)
{
	phy$edge.length<-phy$edge.length^kappa
	return(phy)	
}