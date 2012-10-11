.area.between.curves <-
function(x, f1, f2, xrange=c(0,1))
{
	a<-0.0;
	for(i in 1:length(x)) {
		if(x[i]>=xrange[1] & x[i]<=xrange[2]) {
			if(i==1) {
				lhs<-0  
			} else if(x[i-1]<xrange[1]) {
				lhs<-xrange[1]
			} else lhs<-x[i-1];
			if(i==length(x)) { 
				rhs<-x[i]
			} else if(x[i+1]>xrange[2]) {
				rhs<-xrange[2]; 
			} else rhs<-x[i+1];
			a<-a+(f2[i]-f1[i])*(rhs-lhs)/2;	
		} else if(i!=1) if(x[i-1]>=xrange[1] & x[i-1]<=xrange[2]) {
			y1<-f1[i-1]+(f1[i]-f1[i-1])*(xrange[2]-x[i-1])/(x[i]-x[i-1])
			y2<-f2[i-1]+(f2[i]-f2[i-1])*(xrange[2]-x[i-1])/(x[i]-x[i-1])
			a<-a+(y2-y1)*(xrange[2]-x[i-1])/2;
		} else if(i!=length(x)) if(x[i+1]>=xrange[1] & x[i+1]<=xrange[2]) {
			y1<-f1[i]+(f1[i+1]-f1[i])*(xrange[1]-x[i])/(x[i+1]-x[i])
			y2<-f2[i]+(f2[i+1]-f2[i])*(xrange[1]-x[i])/(x[i+1]-x[i])
			
			a<-a+(y2-y1)*(x[i+1]-xrange[1])/2;
		}
	}
	return(a)
}



disparity <- function(phy=NULL, data, disp=c("avg.sq", "avg.manhattan", "num.states")){
	
	if(!is.null(phy)){
		if (!"phylo"%in%class(phy)) stop("Supply 'phy' as a 'phylo' object")
		td<-treedata(phy, data)
		phy=td$phy
		data=td$data
		desc=.cache_tree(phy)$tips
		nb.tip <- length(td$phy$tip.label)
		nb.node <- td$phy$Nnode
		result<-numeric()
		
		for(i in 1:nb.node) {
			l<-desc[[nb.tip+i]]
			d<-td$data[phy$tip.label[l],]
			result[i]<-.disparity(d, disp)
		}
		names(result)=nb.tip+1:nb.node
		return(result)	  	
		
	} else {
		return(.disparity(data, disp))
	}
	
}


.disparity <- function(data, disp=c("avg.sq", "avg.manhattan", "num.states")){
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

vcv.phylo=function(phy, data=NULL, ...){
	if(is.null(data)) return(.vmat(phy)) else return(.ic.sigma(phy, data))
}

.ic.sigma <-
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


.dtt <-
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
		d<-disparity(phy2, td$data, disp=disp)
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
			d<-disparity(phy2, pp, disp=disp)
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

.dtt.polygon=function(mat, t, alpha=0.05){
	k=nrow(mat)
	dd=c(alpha/2, 1-alpha/2)
	
	tmp=sapply(1:k, function(x) quantile(mat[x,], probs=dd, na.rm=TRUE))
	yy=c(tmp[1,], rev(tmp[2,]))
	xx=c(t, rev(t))
	return(cbind(x=xx, y=yy))
}


dtt<-function(phy, data, disp=c("avg.sq", "avg.manhattan", "num.states"), nsims=1000, mdi.range=c(0,1), plot=TRUE){
	disp=match.arg(disp, c("avg.sq", "avg.manhattan", "num.states"))
	
	td<-treedata(phy, data)
	
	dtt.data<-.dtt(td$phy, td$data, disp=disp)
	ltt<-sort(branching.times(td$phy), decreasing=TRUE)
	ltt<-c(0, (max(ltt)-ltt)/max(ltt));
	
	s<-ic.sigma(td$phy, td$data)
	dtt.sims=NULL
	MDI=NULL
	
	if(is.numeric(nsims)){
		if(nsims>0){
			sims<-sim.char(td$phy, s, nsims)
			
			dtt.sims<-.dtt(td$phy, sims)
			mean.sims<-apply(dtt.sims, 1, mean)
			
			MDI<-unname(.area.between.curves(ltt, apply(dtt.sims, 1, median), dtt.data, mdi.range))
			colnames(dtt.sims)=NULL
		}
	}
	
	if(plot){
#		ylim=c(range(pretty(c(dtt.data, c(dtt.sims)))))
		ylim=c(range(pretty(dtt.data)))
		plot(ltt, dtt.data, xlab="relative time", ylab="disparity", ylim=ylim, bty="n", type="n")
		if(!is.null(dtt.sims)){
#			poly=.dtt.polygon(dtt.sims, ltt, alpha=0.05)
#			polygon(poly[,"x"], poly[,"y"], col=.transparency("lightgray", 0.5), border=NA) 
			lines(ltt, mean.sims, lty=2)
			lines(ltt, dtt.data, type="l", lwd=2)
		}
	}
	
	res=list(dtt=dtt.data, dtt_sims=dtt.sims, times=ltt, MDI=MDI)
	drp=sapply(res, function(x) is.null(x))
	if(any(drp)) res=res[-which(drp)]
	return(res)
}

