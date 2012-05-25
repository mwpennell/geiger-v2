`rate.estimate` <-
function(time=0, n=0, phy=NULL, epsilon = 0, missing = 0, crown=TRUE, kendall.moran=FALSE)
{
	if(!is.null(phy)) {
		if (class(phy) != "phylo") 
       		stop("object \"phy\" is not of class \"phylo\"")
       	
		if(kendall.moran)
		{
			if(missing != 0)
				cat("WARNING: current implementation of Kendall-Moran estimate does not account for missing taxa\n")
			return(kendallmoran.rate(phy))
		}
    	time<-max(branching.times(phy))
		n<-length(phy$tip.label)
	
		n<-n+missing;
		crown=TRUE;
	}
	
    if(crown==TRUE) {
    	if(epsilon==0) {
			rate=(log(n)-log(2))/time
    	} else {
 			rate=1/time*(log(n/2*(1-epsilon^2)+
 					2*epsilon+1/2*(1-epsilon)*sqrt(n*(n*epsilon^2-
 					8*epsilon+2*n*epsilon+n)))-log(2))
    	}
    
    } else {
    	if(epsilon==0) {
			rate=log(n)/time
    	} else {
 			rate=1/time*log(n*(1-epsilon)+epsilon)
    	}
    
    
    }

   return(rate)
}

'kendallmoran.rate' <-
function(phy)
{
	s<-sum(phy$edge.length)
	rate<-(length(phy$tip.label)-2)/s
	return(rate)
}

'crown.p' <- 
function(time, r, epsilon, n)
{
b<-((exp((r*time)))-1)/((exp((r*time)))-epsilon)
a<-epsilon*b
p<-(((b^(n-2)))*((n*(1-a-b+(a*b)))+a+(2*b)-1))/(1+a)
return(p)
}

'stem.p' <-
function(time, r, epsilon, n)
{
b<-((exp((r*time)))-1)/((exp((r*time)))-epsilon)
p<-(b^(n-1))
return(p)
}

'stem.limits' <-
function(r, epsilon, time, prob=c(0.025, 0.975))
{
	limits <- matrix(nrow=length(time), ncol=2)
	for (i in 1:length(time)) {
		beta <- (exp(r*time[i])-1)/(exp(r*time[i]) - epsilon) #From M&S '01 2b
		alpha <- epsilon * beta #From M&S '01 2a
		u <- (log(beta) + log(prob[1]))/log(beta) #From M&S '01 10a
		l <- (log(beta) + log(prob[2]))/log(beta)
		limits[i, 1] <- l
		limits[i, 2] <- u
	}
	return(limits)
}

'crown.limits' <- 
function(r, epsilon, time, prob=c(0.025, 0.975))
{
	limits <- matrix(nrow=length(time), ncol=2)
	for (i in 1:length(time))
	{
		foo<-function(x, prob)
			(crown.p(time[i], r, epsilon, x)-prob)^2
		u<-optim(foo, par=exp(r*time), prob=prob[1], method="L-BFGS-B")
		l<-optim(foo, par=exp(r*time), prob=prob[2], method="L-BFGS-B")

		limits[i, 1] <- l$par
		limits[i, 2] <- u$par
	}
	return(limits)
}

