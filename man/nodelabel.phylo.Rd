\name{nodelabel.phylo}
\alias{nodelabel.phylo}
\alias{lookup.phylo}
\alias{phylo.lookup}
\alias{phylo.clades}
\alias{glomogram.phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Blending information from taxonomies and trees}
\description{
working with systematic reference tables and phylogenies}
\usage{
nodelabel.phylo(phy, taxonomy, strict=TRUE, ncores=NULL)
phylo.lookup(taxonomy, ncores=NULL)
lookup.phylo(phy, taxonomy = NULL, clades = NULL, ncores=NULL)
phylo.clades(clades, phy=NULL, unplaced=TRUE, ncores=NULL)
glomogram.phylo(phy, subtrees)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo' ('multiPhylo' in \code{phylo.clades})}
  \item{taxonomy}{a linkage table (of class \code{matrix}) between tips of the phylogeny and clades represented in the tree; rownames of 'taxonomy' should be tips found in the phylogeny}
  \item{clades}{a named list of clade definitions (i.e., spanning taxa; see \bold{Examples}); spanning taxa may invoke other definitions found within the \code{clades} list}
  \item{unplaced}{whether to use 'unplaced' taxa if given as an element in \code{clades}} 
  \item{subtrees}{a list of trees to be grafted into \code{phy}; each \code{subtrees} element must either be a 'phylo' or 'multiPhylo' object with branch lengths in the same units as \code{phy} (see \bold{Examples})}
  \item{strict}{whether to enforce strict labeling of nodes or allow
    liberal estimates of the best location of a given nodelabel}
  \item{ncores}{the maximum number of cores to be used}
}

\details{
	\code{nodelabel.phylo} provides a function (as part of the \code{phylo} object returned) to resolve the hash key and node identifier for a label found in \code{taxonomy}. 
	This function is the \code{FUN} element of the returned object. If the taxonomic label cannot be properly placed in the tree (i.e., no subtree is found that is absolutely 
	consistent with the supplied \code{taxonomy}, the nearest matching node(s) will be returned when invoking \code{FUN}).   
	
	\code{phylo.lookup} converts a \code{taxonomy} into a phylogenetic tree. 
	
	\code{lookup.phylo} converts a phylogenetic tree (\code{phy}) into a linkage table based on nodelabels associated with \code{phy}, which can be supplemented with a \code{taxonomy} and (or) 
	\code{clades} object.   
	
	\code{phylo.clades} returns a series of phylogenetic subtrees based on clade definitions found in the \code{clades} object.  Definitions can be handles that are recursive (see \bold{Examples}).
}
\author{JM Eastman}

\examples{
\dontrun{
sal=get(data(caudata))
print(head(sal$tax))

## TREE from TABLE: phylo.lookup()
tax=cbind(sal$tax[,c("subfamily", "family", "suborder")], order="Caudata")
tphy=phylo.lookup(tax, ncores=2)
print(tphy)
head(tphy$node.label)

## TABLE from TREE: lookup.phylo()
tax=sal$tax[,c("genus", "family")]
cld=list(
	Sirenoidea=c("Siren", "Pseudobranchus"), 
	Salamandroidea=c("Ambystomatidae", "Plethodontidae"), 
	Cryptobranchoidea=c("Hynobius_naevius", "Cryptobranchus_alleganiensis"),
	CAUDATA=c("Sirenoidea","Salamandroidea","Cryptobranchoidea")
)
lkp=lookup.phylo(sal$phy, taxonomy=tax, clades=cld, ncores=2)
print(lkp)
nphy=nodelabel.phylo(sal$phy, lkp, ncores=2)
dev.new()
plot.phylo(ladderize(nphy,right=FALSE), cex=0.35, 
type="fan", label.offset=2.5, no.margin=TRUE, edge.color="gray", edge.width=0.5)

nodelabels(nphy$node.label, cex=0.45, col="red", frame="n")

## CLADES to TREE: phylo.clades()
fmrca=geiger:::.mrca
salamandroidea=extract.clade(nphy, fmrca(c("Ambystomatidae", "Plethodontidae"), nphy))
cryptobranchoidea=extract.clade(nphy, fmrca(c("Cryptobranchidae", "Hynobiidae"), nphy))
siren=extract.clade(nphy, fmrca(c("Siren_lacertina", "Siren_intermedia"), nphy))

clades=list(
			 Sirenoidea=c("Siren", "Pseudobranchus"), 
			 Caudata=c("Sirenoidea","Salamandroidea","Cryptobranchoidea"),
			 AMPHIBIA=c("Caudata","Anura","Gymnophiona")
)

phy=list(Cryptobranchoidea=cryptobranchoidea, Salamandroidea=salamandroidea, Siren=siren)
class(phy)="multiPhylo"

res=phylo.clades(clades, phy, ncores=2)
amph=nodelabel.phylo(res$AMPHIBIA, lkp, ncores=2)
print(amph$FUN("Salamandroidea"))
dev.new()
plot(ladderize(amph, right=FALSE), cex=0.2, label.offset=0.05)
nodelabels(amph$node.label, cex=0.35, col="red", frame="n")


## GLOMOGRAM
sirenidae=extract.clade(nphy, fmrca(c("Siren_lacertina", "Pseudobranchus_axanthus"), nphy))
ambystomatidae=extract.clade(nphy, fmrca(c("Ambystoma_gracile", "Ambystoma_texanum"), nphy))
trees=list(
	Cryptobranchoidea=cryptobranchoidea,
	Sirenidae=sirenidae,
	Ambystomatidae=ambystomatidae
)
class(trees)="multiPhylo"

fam=sal$fam
ftax=unique(sal$tax[,c("family", "suborder")])
rownames(ftax)=unname(ftax[,"family"])
fam=nodelabel.phylo(fam, ftax, ncores=2)
fam$FUN("Salamandroidea")

res=glomogram.phylo(fam, trees)
dev.new()
zz=match(res$tip.label, fam$tip.label)
cc=integer(length(zz))
cc[!is.na(zz)]=1

plot(ladderize(res, right=FALSE), cex=1, label.offset=5, tip.color=ifelse(cc==1, "red", "black"))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data}
\keyword{graphs}
\keyword{manip}

