`get.simulation.matrix` <-
function(phy)
{
	if (class(phy) != "phylo")
		stop("object \"phy\" is not of class \"phylo\"");
	tmp<-as.numeric(phy$edge)
	nb.tip <- max(tmp)
    bl <- phy$edge.length;
    result<-matrix(0, nrow=nb.tip, ncol=length(bl));

	for(i in 1:nb.tip) {
		j<-i;
		while(j != "-1") {
			branch<-which(as.numeric(phy$edge[,2]) == j);
			result[i, branch]<-sqrt(bl[branch]);
			j<-phy$edge[branch, 1];
		}
	}
	
	return(as.matrix(result));
}

