prune.extinct.taxa<-function (phy, tol = .Machine$double.eps^0.5) 
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

prune.random.taxa<-function (phy, n) 
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