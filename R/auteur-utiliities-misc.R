print.rbm<-function (x, ...) 
{
    cat("\nBrownian Motion likelihood function (1D relaxed):\n\n")
	attributes(x)=NULL
    print(unclass(x))
}

print.rbm=function (x, printlen = 3, ...) 
{
    cat("likelihood function for relaxed-rates branching diffusion\n")
	aa=names(argnames(x))
    cat("\targnames:", paste(aa, collapse = ", "))		 
	cat("\n\n")
	f=x
	attributes(f)=NULL
	cat("definition:\n")
	print(f)
	
	cat("\n\nNOTE: 'rates' vector should be given in order given by the function 'argnames()'")
}


.check.multicore=function(){
	tmp=rownames(installed.packages())
	if("multicore"%in%tmp) {
		require(multicore)
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

