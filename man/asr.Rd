\name{asr}
\alias{asr}
\title{estimate ancestral character states for continuous characters}
\description{
	estimating ancestral character states for continuous characters
}
\usage{
asr(phy, dat) 
}
\arguments{
  \item{dat}{ data vector for a single continuous trait, with names matching tips in \code{phy}}
  \item{phy}{a phylogenetic tree of class 'phylo'}
}
\details{
	This function uses an algorithm to calculate Maximum-likelihood estimates of 
	ancestral character states for continuous characters.
	
	The \pkg{ape} function \code{\link[ape]{ace}} carries out a similar function but finds estimates by actually maximizing the likelihood across all ancestors simultaneously. The algorithm used here is faster and sometimes more reliable especially for larger trees.
    The \pkg{phytools} function \code{\link[phytools]{fastAnc}} carries out this function
}

\value{
	The function returns a vector of ancestral character states at each node, named based on the numeric 
	node-identifiers used inside the \code{phy$edge} object.
	
	If multiple characters are given as \code{dat}, a matrix of ancestral states is returned (note that each 
	character is assumed independent of all others).  
}
	
\author{Luke J. Harmon}

\seealso{
\code{\link[ape]{ace}} in the \pkg{ape} package
}

\examples{

tmp=get(data(geospiza))
geo=treedata(tmp$phy, tmp$dat)

asr(geo$phy, geo$data[,1])

asr(geo$phy, geo$data)
 

}
\keyword{arith}
