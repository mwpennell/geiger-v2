`node.leaves` <-
function(phy, node)
{
	node<-as.numeric(node)
	n<-length(phy$tip.label);
	if(node <= n) return(phy$tip.label[as.numeric(node)])

	l<-character();
	
	d<-node.sons(phy, node);
	for(j in d) {
		if(j <= n) l<-c(l, phy$tip.label[as.numeric(j)])
		else l<-c(l, node.leaves(phy, j));
	}
	return(l);
}

'node.sons' <-
function(phy, node)
{
	r<-which(phy$edge[,1]==node)
	return(phy$edge[r,2])	
}