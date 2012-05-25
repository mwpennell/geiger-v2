`tip.disparity` <-
function(phy, data, data.names=NULL, disp="avg.sq")
{
	if (class(phy) != "phylo")
		stop("object \"phy\" is not of class \"phylo\"");
	
	td<-treedata(phy, data, data.names)

		
	nb.tip <- length(td$phy$tip.label)
    nb.node <- td$phy$Nnode
    result<-numeric();
    
    for(i in 1:nb.node) {
    	l<-node.leaves(td$phy, nb.tip+i);
    	d<-td$data[match(l, row.names(data)),];
    	result[i]<-disp.calc(d, disp);
    }
    return(result);	  	
}

