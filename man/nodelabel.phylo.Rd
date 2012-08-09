\name{nodelabel.phylo}
\alias{nodelabel.phylo}
\alias{lookup.phylo}
\alias{phylo.lookup}
\alias{phylo.clades}
\alias{bind.phylo}
\alias{glomogram.phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Blending information from taxonomies and trees}
\description{
Working with systematic reference tables and phylogenies}
\usage{
nodelabel.phylo(phy, taxonomy, strict=TRUE)
phylo.lookup(taxonomy)
lookup.phylo(phy, taxonomy = NULL, clades = NULL)
phylo.clades(clades, phy=NULL, unplaced=TRUE)
bind.phylo(phy, taxonomy)
glomogram.phylo(phy, subtrees)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo' ('multiPhylo' in \code{phylo.clades})}
  \item{taxonomy}{a ('matrix') linkage table between tips of the phylogeny and clades represented in the tree; rownames of 'taxonomy' should be tips found in the phylogeny}
  \item{clades}{a named list of clade definitions (i.e., spanning taxa; see \bold{Examples}); spanning taxa may invoke other definitions found within the \code{clades} list}
  \item{unplaced}{whether to use 'unplaced' taxa if given as an element in \code{clades}} 
  \item{subtrees}{a list of trees to be grafted into \code{phy}; each \code{subtrees} element has two sub-elements: \code{subtree} and \code{age} (see \bold{Examples})}
  \item{strict}{whether to enforce strict labeling of nodes or allow liberal estimates of the best location of a given node label}
}

\details{
	\code{nodelabel.phylo} provides a function (as part of the \code{phylo} object returned) to resolve the hash key and node identifier for a label found in \code{taxonomy}. 
	This function is the \code{FUN} element of the returned object. If the taxonomic label cannot be properly placed in the tree (i.e., no subtree is found that is absolutely 
	consistent with the supplied \code{taxonomy}, the nearest matching node(s) will be returned when invoking \code{FUN}).     
}
\author{JM Eastman}

\examples{
sal=get(data(urodela))
print(head(sal$tax))

## TREE from TABLE: phylo.lookup()
tax=cbind(sal$tax[,c("subfamily", "family", "suborder")], order="Caudata")
tphy=phylo.lookup(tax)
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
lkp=lookup.phylo(sal$phy, taxonomy=tax, clades=cld)
print(lkp)
nphy=nodelabel.phylo(sal$phy, lkp)
dev.new()
plot.phylo(ladderize(nphy,right=FALSE), cex=0.35, type="fan", label.offset=2.5, no.margin=TRUE, edge.color="gray", edge.width=0.5)

nodelabels(nphy$node.label, cex=0.45, col="red", frame="n")

## BIND TAXA to TREE: bind.phylo()
bphy=bind.phylo(sal$fam, sal$tax)
bphy=nodelabel.phylo(bphy, sal$tax)
dev.new()
plot.phylo(ladderize(bphy,right=FALSE), cex=0.35, show.tip=FALSE, label.offset=2.5, no.margin=TRUE, edge.color="gray", edge.width=0.5)

nodelabels(bphy$node.label, cex=0.45, col="red", frame="n")

## CLADES to TREE: phylo.clades()
salamandroidea=extract.clade(nphy, geiger:::.mrca(c("Ambystomatidae", "Plethodontidae"), nphy))
cryptobranchoidea=extract.clade(nphy, geiger:::.mrca(c("Cryptobranchidae", "Hynobiidae"), nphy))
siren=extract.clade(nphy, geiger:::.mrca(c("Siren_lacertina", "Siren_intermedia"), nphy))

clades=list(
			 Sirenoidea=c("Siren", "Pseudobranchus"), 
			 Caudata=c("Sirenoidea","Salamandroidea","Cryptobranchoidea"),
			 AMPHIBIA=c("Caudata","Anura","Gymnophiona")
)

phy=list(Cryptobranchoidea=cryptobranchoidea, Salamandroidea=salamandroidea, Siren=siren)
class(phy)="multiPhylo"

res=phylo.clades(clades, phy)
amph=nodelabel.phylo(res$AMPHIBIA, lkp)
dev.new()
plot(ladderize(amph, right=FALSE), cex=0.2, label.offset=0.05)
nodelabels(amph$node.label, cex=0.35, col="red", frame="n")


## GLOMOGRAM
sirenidae=extract.clade(nphy, geiger:::.mrca(c("Siren_lacertina", "Pseudobranchus_axanthus"), nphy))
ambystomatidae=extract.clade(nphy, geiger:::.mrca(c("Ambystoma_gracile", "Ambystoma_texanum"), nphy))
trees=list(
	Cryptobranchoidea=cryptobranchoidea,
	Sirenidae=sirenidae,
	Ambystomatidae=ambystomatidae
)
class(trees)="multiPhylo"

fam=sal$fam
ftax=unique(sal$tax[,c("family", "suborder")])
rownames(ftax)=unname(ftax[,"family"])
fam=nodelabel.phylo(fam, ftax)
fam$FUN("Salamandroidea")

res=glomogram.phylo(fam, trees)
dev.new()
zz=match(res$tip.label, fam$tip.label)
cc=integer(length(zz))
cc[!is.na(zz)]=1

plot(ladderize(res, right=FALSE), cex=1, label.offset=5, tip.color=ifelse(cc==1, "red", "black"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
