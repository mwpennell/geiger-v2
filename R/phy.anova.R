phy.anova<-function(phy, data, group, data.names=NULL, nsim=1000)
{
	td<-treedata(phy, data, data.names)
	
	s<-mean(pic(td$data, td$phy)^2)
	a<-anova(lm(td$data~group))
	f.data<-a[1,4]
	sims<-sim.char(td$phy, as.matrix(s), nsim=nsim)
	
	foo<-function(xx) anova(lm(xx~group))[1,4]
	f.null<-apply(sims, 3, foo)

	cat("Standard ANOVA:\n")	
	print(a)

	cat("\n\nPhylogenetic p-value: \t")
	cat((sum(f.null>f.data)+1)/(nsim+1))

		
}

phy.manova<-function(phy, data, group, data.names=NULL, nsim=1000, test="Wilks")
{	
	td<-treedata(phy, data, data.names)

	s<-ic.sigma(td$phy, td$data)
	
	m<-summary.manova(manova(as.matrix(td$data)~group), test=test)
	
	w.data<-m[[4]][1,2]
	
	sims<-sim.char(td$phy, s, nsim=nsim)
	
	foo<-function(xx) summary.manova(manova(as.matrix(xx)~group), test=test)[[4]][1,2]
	
	w.null<-apply(sims, 3, foo)

	cat("Standard MANOVA:\n")	
	print(m)
	
	cat("\n\nPhylogenetic p-value: \t")
	cat((sum(w.data>w.null)+1)/(nsim+1))
	
}