\name{dcount}
\alias{dcount}
\title{
prior densities for truncated discrete random variable
}
\description{
creating a prior density function for a truncated discrete random variable}
\usage{
dcount(x, FUN, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an integer vector spanning from the minimal to maximal values}
  \item{FUN}{a probability density function (see \bold{Examples})}
  \item{\dots}{additional arguments to be passed to \code{FUN}}
}

\seealso{\code{\link[geiger]{make.gbm}}}

\author{JM Eastman}

\examples{
range=0:100
u=dcount(range, FUN=dunif, min=min(range), max=max(range))
g=dcount(range, FUN=dgeom, prob=0.05)
p=dcount(range, FUN=dtpois, min=min(range), max=max(range), lambda=0.5)
priors=list(pois=p, geom=g, unif=u)
plot(range(range), c(0,1), xlab="index", ylab="cumsum(prob)", type="n", bty="n")
for(i in 1:length(priors)){
	points(attributes(attributes(priors[[i]])$density)$cumsum, col=i, pch=22, cex=0.35)
}
legend("bottomright", bty="n", legend=names(priors), col=1:length(priors), pch=22)


## LN prior probabilities
print(u(0))	## dunif
print(g(0))	## dgeom
print(p(0)) ## dtpois
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arith}
