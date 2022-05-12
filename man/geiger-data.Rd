\name{geiger-data}
\alias{amphibia}
\alias{chelonia}
\alias{geospiza}
\alias{primates}
\alias{caudata}
\alias{whales}
\alias{carnivores}
\alias{caniformia}

\title{example datasets}
\description{providing access to comparative datasets}
\usage{
data(amphibia)
data(caniformia)
data(carnivores)
data(caudata)
data(chelonia)
data(geospiza)
data(primates)
data(whales)
}

\details{
The objects \code{caudata}, \code{chelonia}, \code{carnivores}, \code{geospiza}, and \code{primates} each have at least two items, a \code{phy} object 
and a \code{dat} object. The \code{phy} object is a phylogenetic tree of class 'phylo' (see \code{\link[ape]{read.tree}}).  The \code{dat} object (e.g., \code{caudata$dat})
is a named vector of (natural log-transformed) body sizes for each group. The salamander data (object \code{caudata}) also includes 
the systematics for all recognized taxa (object \code{caudata$tax}) as well as a time-calibrated family-level phylogeny (\code{caudata$fam}) from Zhang and Wake (2009). The object \code{caudata$phy} is an unpublished glomogram of mostly family level phylogenies from the literature.  The backbone for that tree is from Zhang and Wake (2009). 

The \code{amphibia} object is a set of three trees, the last of which is a time-scaled estimate of the Pyron and Wiens (2011) tree (see \code{\link[geiger]{congruify.phylo}}).

The \code{whales} object is a dataset including a tree and a taxon richness matrix (see \code{\link[geiger]{medusa}}).
}

\references{

Data are from the following sources:

\code{AMPHIBIANS (amphibia)}

\itemize{
		\item{Roelants K, DJ Gower, M Wilkinson, SP Loader, SD Biju, K Guillaume, L Moriau, and F Bossuyt. 2007. Global patterns of diversification in the history of modern amphibians. \emph{PNAS} 104:887-892.}
		
		\item{Pyron RA and JJ Wiens. 2011. A large-scale phylogeny of Amphibia including over 2800 species, and a revised classification of extant frogs, salamanders, and caecilians. \emph{MPE} 61:543-583.}
}

\code{SALAMANDERS (caudata)}

\itemize{

		\item{Adams DC, CM Berns, KH Kozak, and JJ Wiens. 2009. Are rates of species diversification correlated with rates of morphological evolution? \emph{PRSB} 276:2729-2738. }
		
		\item{Bonett RM, PT Chippindale, PE Moler, RW van Devender, and DB Wake. 2009. Evolution of gigantism in amphiumid salamanders. \emph{PLoSONE} 4(5):e5615.}
				
		\item{Kozak KH, RW Mendyk, and JJ Wiens. 2009. Can Parallel Diversification Occur in Sympatry? Repeated Patterns of Body-Size Evolution in Coexisting Clades of North American Salamanders. \emph{Evolution} 63:1769-1784.}
		
		\item{Weisrock DW, TJ Papenfuss, JR Macey, SN Litvinchuk, R Polymeni, IH Ugurtas, E Zhao, H Jowkar, and A Larson. 2006. A molecular assessment of phylogenetic relationships and lineage accumulation rates within the family Salamandridae (Amphibia, Caudata). \emph{MPE} 41:368-383.}
		
		\item{Wiens JJ and JT Hoverman. 2008. Digit reduction, body size, and paedomorphosis in salamanders. \emph{Evolution and Development} 10:449-463.}
		
		\item{Zhang P, Y-Q Chen, H Zhou, X-L Wang, TJ Papenfuss, DB Wake and L-H Qu. 2006. Phylogeny, evolution, and biogeography of Asiatic salamanders (Hynobiidae). \emph{PNAS} 103:7360-7365.}
		
		\item{Zhang P and DB Wake. 2009. Higher-level salamander relationships and divergence dates inferred from complete mitochondrial genomes. \emph{MPE} 53:492-508. }
}

\code{PRIMATES (primates)}

\itemize{
		\item{Redding DW, C DeWolff, and AO Mooers. 2010. Evolutionary distinctiveness, threat status and ecological oddity in primates. \emph{Conservation Biology} 24:1052-1058.}
		
		\item{Vos RA and AO Mooers. 2006. A new dated supertree of the Primates. Chapter 5. In: VOS RA (Ed.) \emph{Inferring large phylogenies: the big tree problem}. [Ph.D. thesis]. Burnaby BC, Canada: Simon Fraser University.}
}
	
\code{CARNIVORES (carnivores)}

\itemize{
		\item{Eizirik E, WJ Murphy, K-P Koepfli, WE Johnson, JW Dragoo, RK Wayne, and SJ O'Brien. 2010. Pattern and timing of diversification of the mammalian order Carnivora inferred from multiple nuclear gene sequences. \emph{Molecular Phylogenetic and Evolution} 56:49-63.}
        
        \item{Wozencraft WC. 2005. Order Carnivora in Wilson DE and DM Reeder (Eds.) \emph{Mammal Species of the World}. Johns Hopkins University Press.}
         
        \item{Jones KE, J Bielby, M Cardillo, et al. 2009. \emph{Ecological Archives} E090-184.}
}	


\code{CANIFORMS (caniformia)}

\itemize{

		\item{Slater GJ, LJ Harmon, and ME Alfaro. 2012. Integrating fossils with molecular phylogenies improves inference of trait evolution. \emph{Evolution} 66:3931-3944. }
}
	
\code{TURTLES (chelonia)}

\itemize{
		\item{Jaffe AL, GJ Slater, and ME Alfaro. 2011. Ecological habitat and body size evolution in turtles. \emph{Biology Letters} 7:558-561.}
}	


\code{DARWINS FINCHES (geospiza)}

\itemize{

		\item{Schluter D}

}

\code{WHALES (whales)}

\itemize{
		\item{data compiled by GJ Slater}
		
		\item{Paleobiology Database. 2011. \url{https://paleobiodb.org/}.}
}

}
\keyword{datasets}

