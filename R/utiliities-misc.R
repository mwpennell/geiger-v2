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



#general printing utility for ensuring equal numbers of characters within columns and defining spacing between columns
#author: JM EASTMAN 2010
#note: works only for numeric dataframes

.print.table=function(df,digits=4,buffer=5){
	if(length(buffer) != ncol(df) | length(buffer)==1) buffer=rep(buffer[1],ncol(df))
	if(length(digits) != ncol(df) | length(digits)==1) digits=rep(digits[1],ncol(df))
	ss=sapply(round(df),nchar)
	lar=df>1
	nn=sapply(names(df),nchar)
	
    # find longest string
	strw=sapply(1:ncol(df), function(x) max(nn, max(1,(ss[lar])+digits[x],na.rm=TRUE),na.rm=TRUE))
	pr.df=data.frame(do.call(cbind, lapply(1:ncol(df), function(x) sprintf(paste("%",(strw[x]+buffer[x]),".",digits[x],"f",sep=""),df[,x]))))
	names(pr.df)=names(df)
	rownames(pr.df)=rownames(df)
	print(pr.df)
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

