coef.gfit=function(object, ...){
    if(is.constrained(object$lik)) p=names(object$lik(argn(object$lik),pars.only=TRUE)) else p=argn(object$lik)
    unlist(object$opt[p])
}

coef.gfits=function(object, ...){
    lapply(object, coef)
}

logLik.gfit=function(object, ...){
    object$opt$lnL
}

logLik.gfits=function(object, ...){
    lapply(object, function(x) x$opt$lnL)
}




.get.parallel=function(...){
	
	if(.check.parallel()) {
		if(Sys.getenv("R_PARALLEL")=="FALSE") {
            fx=function(X,FUN,...) lapply(X,FUN,...)
        } else {
            fx=function(X,FUN,...) mclapply(X,FUN,...,mc.silent=TRUE)
        }
	} else {
        fx=function(X,FUN,...) lapply(X,FUN,...)
	}
	
	fx
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

