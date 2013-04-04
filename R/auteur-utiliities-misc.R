print.rbm=function (x, printlen = 3, ...) 
{
    cat("likelihood function for relaxed-rates univariate continuous trait evolution\n")
	aa=names(argn(x))
    cat("\targument names:", paste(aa, collapse = ", "))		 
	cat("\n\n")
	f=x
	attributes(f)=NULL
	cat("definition:\n")
	print(f)
	
	cat("\n\nNOTE: 'rates' vector should be given in order given by the function 'argn()'")
}


.get.parallel=function(){
	ee=Sys.getenv()
	
	if(.check.parallel()) {
		if("ignoreMULTICORE"%in%names(ee)) f=lapply else f=function(X,FUN) mclapply(X,FUN,mc.silent=TRUE)
	} else {
		f=lapply
	}
	
	f
}

.check.parallel=function(){
	tmp=rownames(installed.packages())
	if("parallel"%in%tmp) {
		require(parallel)
		return(TRUE)
	} else {
		return(FALSE)
	}
}




.transparency <- function (col, alpha) 
{
    tmp <- col2rgb(col)/255
    rgb(tmp[1, ], tmp[2, ], tmp[3, ], alpha = alpha)
}

.withinrange <- function (x, min, max) 
{
    a = sign(x - min)
    b = sign(x - max)
    if (abs(a + b) == 2) 
	return(FALSE)
    else return(TRUE)
}

.basename.noext=function(path=""){
	return(sub("[.][^.]*$", "", basename(path), perl=TRUE))
	
}



.resolve.executable=function(package="geiger"){
	packagedir=system.file(package=package)
	execs=lapply(d<-dir(paste(packagedir,"exec",sep="/")), function(x) {paste(packagedir, "exec", x, sep="/")})
	names(execs)=.basename.noext(d)
	return(execs)
}

