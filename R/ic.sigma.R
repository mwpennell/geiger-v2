`ic.sigma` <-
function(phy, data, data.names=NULL)
{
	td<-treedata(phy, data, data.names, sort=T)
	
	f<-function(x) pic(x, td$phy)
	ic<-apply(td$data, 2, f)
	r<-crossprod(ic, ic)/nrow(ic)
	return(r)
}

