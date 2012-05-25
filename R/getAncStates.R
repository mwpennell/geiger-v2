getAncStates<-function(x, phy) {
	
	nn<-phy$Nnode
	nt<-length(phy$tip.label)
	anc<-numeric(nn)
	
	for(i in 1:nn) {
		node<-nt+i
		rp<-root(phy, node=node)
		rp<-multi2di(rp)
		pp<-picAnc(x, rp)
		mm<-match(rp$edge[,1], rp$edge[,2])
		rr<-which(is.na(mm))[1]
		root<-rp$edge[rr,1]
		anc[i]<-pp[root]
		}
		
	names(anc)<-1:nn + nt	
	anc
	
	}



picAnc<-function(x,phy,sd=NULL,n=NULL,se=NULL,scaled=TRUE,var.contrasts=FALSE)
{
	if (class(phy) != "phylo") 
        stop("object 'phy' is not of class \"phylo\"")
    if (is.null(phy$edge.length)) 
        stop("your tree has no branch lengths: you may consider setting them equal to one, or using the function `compute.brlen'.")
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
   
    if (length(x) != nb.tip) 
        stop("length of phenotypic and of phylogenetic data do not match")
    if (any(is.na(x))) 
        stop("the present method cannot (yet) be used directly with missing data: you may consider removing the species with missing data from your tree with the function `drop.tip'.")
    phy<-reorder(phy, "pruningwise")

    phenotype<-ssq<-deltaV<-numeric(nb.tip + nb.node)
  	deltaV[1:nb.tip]<-0
    if (is.null(names(x))){
        phenotype[1:nb.tip] <- x
		if(!is.null(n))ssq[1:nb.tip]<-1/n
	}else{
		if (all(names(x) %in% phy$tip.label)){
			phenotype[1:nb.tip] <- x[phy$tip.label]
			if(!is.null(n))ssq[1:nb.tip]<-1/n[phy$tip.label]
		}else{
			phenotype[1:nb.tip] <- x
            if(!is.null(n))ssq[1:nb.tip]<-1/n
            warning("the names of argument \"x\" and the names of the tip labels did not match: the former were ignored in the analysis.")
        }
    }
    contr<-var.contr<-numeric(nb.node)
#pic function brought in from C source
for(i in seq(1,nb.tip*2-2,by=2)){
	j<-i+1
	anc<-phy$edge[i,1]
	d1<-phy$edge[i,2]
	d2<-phy$edge[j,2]
	sumbl<-phy$edge.length[i]+phy$edge.length[j]
	ic<-anc-nb.tip
	contr[ic]<-phenotype[d1]-phenotype[d2]
	if(scaled)contr[ic]<-contr[ic]/sqrt(sumbl)
	if(var.contrasts)var.contr[ic]<-sumbl
	phenotype[anc]<-(phenotype[d1]*phy$edge.length[j]+phenotype[d2]*phy$edge.length[i])/sumbl
	if(j!=(nb.tip*2-2)){
		k<-which(phy$edge[,2]==anc)
		phy$edge.length[k]<-phy$edge.length[k]+phy$edge.length[i]*phy$edge.length[j]/sumbl
	}
	}
	phenotype
}
