######################################################################################
##STOCHASTIC SIMULATION OF A TIME-HOMOGENOUS BIRTH-DEATH PROCESS
##PARAMETERS: N0=starting number of lineages, ks=speciation rate, ep=relative rate of extinction, finaltime=endpoint for simulation
######################################################################################
BDsim <- function(nStart, b, d, times) {
	n <- nStart
	t <- 0
	pop<-numeric(length(times))
	pop[]=0
	pop[1] <- nStart
	
	times<-c(0, sort(times))
	i=2
	while(1) {
		if(t>times[i]) {
			m=max(which(t>times))
			pop[i:m]=n;
			i=m+1;	
		}

		waittime <- rexp(1, rate=n*(b+d))
		t <- t+ waittime
		if(t > max(times)) {
			break 
		} else {
			ran <- runif(1)
			if(ran<=b/(b+d)) n <- n+1						else n <- n-1
		}
		if(n==0) break
	}
	pop[i:length(times)]=n
	res=cbind(times, pop)
	colnames(res)=c("Time", "n")		
	return(res)
}	

