`disp.calc` <-
function(data, disp=c("avg.sq", "avg.manhattan", "num.states")){
	disp=match.arg(disp, c("avg.sq", "avg.manhattan", "num.states"))
	
	if(disp=="avg.sq") {
		d<-dist(data, method="euclidean")^2
		r<-mean(d)
	}
	else if(disp=="avg.manhattan") {
		d<-dist(data, method="manhattan")
		r<-mean(d)
	}
	else if(disp=="num.states") {
		f<-function(x) length(unique(x))
		d<-apply(data, 2, f)
		r<-mean(d)
	}
	else r<-0;
	return(r)
}

`ic.sigma` <-
function(phy, data)
{
	td<-treedata(phy, data, sort=TRUE)
	
	f<-function(x) pic(x, td$phy)
	ic<-apply(td$data, 2, function(x) {
				names(x)=rownames(td$data)
				f(x)
	})
	r<-crossprod(ic, ic)/nrow(ic)
	return(r)
}



`tip.disparity` <-
function(phy, data,  disp=c("avg.sq", "avg.manhattan", "num.states")){
	disp=match.arg(disp, c("avg.sq", "avg.manhattan", "num.states"))
	
	if (class(phy) != "phylo")
	stop("object \"phy\" is not of class \"phylo\"");
	
	td<-treedata(phy, data)
	
	
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


`dtt` <-
function(phy, data, disp=c("avg.sq", "avg.manhattan", "num.states")){
	disp=match.arg(disp, c("avg.sq", "avg.manhattan", "num.states"))

	phy$node.label<-NULL
	td<-treedata(phy, data)
	phy2<-td$phy
	phy<-new2old.phylo(td$phy)
	
	result<-numeric()
	

	node.depth<-branching.times(phy2);
	stem.depth<-numeric();
	stem.depth[1]<-node.depth[1];
	for(i in 2:phy2$Nnode) {
		anc<-which(as.numeric(phy$edge[,2])==-i)
		stem.depth[i]<-node.depth[names(node.depth)==phy2$edge[anc,1]]
	}
		
	ltt<-sort(node.depth, decreasing=TRUE)
	node.depth<-node.depth/max(ltt);
	stem.depth<-stem.depth/max(ltt);
	ltt<-ltt/max(ltt);
	if(length(dim(td$data))==2) {
		d<-tip.disparity(phy2, td$data, disp=disp);
		result[1]<-d[1]
		for(i in 2:length(ltt)) {
			x<-d[stem.depth>=ltt[i-1]&node.depth<ltt[i-1]]
			if(length(x)==0) result[i]=0
			else result[i]<-mean(x);
		}
		result[length(ltt)+1]<-0;
		if(result[1]>0)
			result<-result/result[1];
			
	} else {
		if(length(dim(td$data))!=3)
			stop("Error in data");
		
		for(i in 1:dim(td$data)[3]) {
			pp<-as.matrix(td$data[,,i])
			d<-tip.disparity(phy2, pp, disp=disp);
			y<-numeric()
	
			y[1]<-d[1]
			for(j in 2:length(ltt)) {
				x<-d[stem.depth>=ltt[j-1]&node.depth<ltt[j-1]]
				if(length(x)==0) y[j]=0
				else y[j]<-mean(x);
			}
			y[length(ltt)+1]<-0;
			if(y[1]>0)
			y<-y/y[1];
			
			result<-cbind(result, y)
		}
	}
	
	return(result);	
}

dtt.full<-function(phy, data, disp=c("avg.sq", "avg.manhattan", "num.states"), nsims=1000, mdi.range=c(0,1)){
	disp=match.arg(disp, c("avg.sq", "avg.manhattan", "num.states"))

	td<-treedata(phy, data)

	dtt.data<-dtt(td$phy, td$data, disp=disp)
	ltt<-sort(branching.times(td$phy), decreasing=TRUE)
	ltt<-c(0, (max(ltt)-ltt)/max(ltt));
	plot(ltt, dtt.data, type="l", lwd=2, xlab="Relative time", ylab="Disparity");
	
	s<-ic.sigma(td$phy, td$data)
	sims<-sim.char(td$phy, s, nsims)

	dtt.sims<-dtt(td$phy, sims)
	mean.sims<-apply(dtt.sims, 1, mean)

	lines(ltt, mean.sims, lty=2)

	MDI<-area.between.curves(ltt, apply(dtt.sims, 1, median), dtt.data, mdi.range)

	return(list(dtt.data=dtt.data, dtt.sims=dtt.sims, times=ltt, MDI=MDI))
	
}

