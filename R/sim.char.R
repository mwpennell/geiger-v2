`sim.char` <-
function(phy, model.matrix, nsims=1, model="brownian", root.state=1)
{
	phy<-new2old.phylo(phy)
	
	nchar<-nrow(model.matrix);
	m<-get.simulation.matrix(phy);
	nbranches<-ncol(m);	
	nspecies<-nrow(m);

	if(model=="brownian" | model=="speciational")
	{
		if(model=="speciational") {
			m[m>0]<-1.0;
		}
	
		rnd<-t(mvrnorm(nsims*nbranches, mu=rep(0, nchar), Sigma=model.matrix));
		rnd<-array(rnd, dim=c(nchar, nbranches, nsims));
		
		simulate<-function(v) m %*% as.matrix(v);
		
		result<-apply(rnd, 1, simulate)
		result<-aperm(array(result, dim=c(nspecies, nsims, nchar)), c(1, 3, 2))
		
		rownames(result)<-phy$tip.label;
	}
	
	if(model=="discrete")
	{
		nchar<-length(model.matrix);
		node.value<-numeric(nbranches)
		result<-array(0, dim=c(nspecies, nchar, nsims))
		for(j in 1:nchar) {
			m<-model.matrix[[j]];
			for(k in 1:nsims) {
	   			for(i in 1:nbranches) {
					if(as.numeric(phy$edge[i,1])==-1) s<-root.state
					else {
						parent<-which(phy$edge[,2]==phy$edge[i,1])
						s<-node.value[parent]
					}	
					p<-MatrixExp(m*phy$edge.length[i])
					probs<-cumsum(p[s,])
					r<-runif(1)
					node.value[i]<-min(which(r<probs))
	   			}
	   			result[,j,k]<-node.value[as.numeric(phy$edge[,2])>0]
	   		}	
		}
		rownames(result)<-phy$tip.label;

	}		
	return(result);
}

