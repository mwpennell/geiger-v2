`set.seed.clock` <-
function(print=F){
	date = date()
 	seed1 = as.numeric(strsplit(substring(date,12,19),":")[[1]])%*%c(1,100,10000)
 	seed <- runif(1, min=0, max=50) * seed1
 	set.seed(seed)
 	if(print) cat("Seed = ", seed, "\n");
 	seed[1,1]
}

