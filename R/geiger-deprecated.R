.deprecate=function(prev, curr){
	warning(paste(sQuote(prev), "is being deprecated: use", paste(sQuote(curr), collapse=" or "), sep=" "))
}


name.check=function(phy, data){
    if(is.null(names(data))) stop("'data' must be given as a vector or matrix with names")
	.deprecate("name.check", "geiger:::.treedata")
	.treedata(phy, data)
}

BDsim=function(nStart, b, d, times){
	.deprecate("BDsim", "sim.bd")
	sim.bd(b=b, d=d, n0=nStart, times=times)
}


birthdeath.tree=function(b, d, time.stop=0, taxa.stop=0, seed=0, print.seed=FALSE, return.all.extinct=TRUE){
	.deprecate("birthdeath.tree", "sim.bdtree")
	crit=.check.stoppingcrit(time.stop, taxa.stop)
	sim.bdtree(b=b, d=d, stop=crit, n=taxa.stop, t=time.stop, seed=seed, extinct=return.all.extinct) 
}

tip.disparity=function(phy, data, disp=c("avg.sq", "avg.manhattan", 
"num.states")){
	.deprecate("tip.disparity", "disparity")
	disparity(phy=phy, data=data, disp=disp)
}

ic.sigma=function(phy, data){
	.deprecate("ic.sigma", "vcv")
	vcv(phy=phy, data=data)
}

rate.estimate=function(time=0, n=0, phy=NULL, epsilon = 0, missing = 0, crown=TRUE, kendall.moran=FALSE){
	.deprecate("rate.estimate", c("bd.ms", "bd.km"))
	if(kendall.moran){
		bd.km(phy=phy, time=time, n=n, missing=missing, crown=crown)
	} else {
		bd.ms(phy=phy, time=time, n=n, missing=missing, crown=crown, epsilon=epsilon)
	}
}

node.leaves=function(phy, node){
	.deprecate("node.leaves", "tips")
	tips(phy, node)

}

getAncStates=function(x, phy){
    .deprecate("getAncStates", "asr")
    asr(phy, x)
}