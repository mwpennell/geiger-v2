\name{auteur.data}
\alias{primates}
\alias{chelonia}
\alias{urodela}
\title{example datasets for auteur}
\description{Provides example body-size data and a corresponding phylogeny}
\usage{data(primates); data(chelonia); data(urodela)}
\details{Each object is a list of at least two items, a \code{phy} object and a \code{dat} object.  The \code{phy} object is a phylogenetic tree of class 'phylo' (see \code{\link[ape]{read.tree}}).  The \code{dat} object 
is a named vector of (natural log-transformed) body sizes for each group. The salamander data (\code{urodela}) also includes the systematics for all recognized taxa as well as a time-calibrated family-level phylogeny from Zhang and Wake (2009).}

\references{ Data are from the following sources:

\code{PRIMATES}

\itemize{
		\item{REDDING DW, C DEWOLFF, and AO MOOERS. 2010. Evolutionary distinctiveness, threat status and ecological oddity in primates. Conservation Biology 24:1052-1058.}
		
		\item{VOS RA and AO MOOERS. 2006. A new dated supertree of the Primates. Chapter 5. In: VOS RA (Ed.) Inferring large phylogenies: the big tree problem. [Ph.D. thesis]. Burnaby BC, Canada: Simon Fraser University.}
}
	
\code{TURTLES}

\itemize{
		\item{JAFFE, G SLATER, and M ALFARO. 2011. Ecological habitat and body size evolution in turtles. Biology Letters: in press. }
}	

\code{SALAMANDERS}

\itemize{

		\item{ADAMS DC, CM BERNS, KH KOZAK, and JJ WIENS. 2009. Are rates of species diversification correlated with rates of morphological evolution? Proceedings of the Royal Society of London B 276:2729-2738. }
		
		\item{BONETT RM, PT CHIPPINDALE, PE MOLER, RW VAN DEVENDER, and DB WAKE. 2009. Evolution of gigantism in amphiumid salamanders. PLoSONE 4(5):e5615.}
				
		\item{KOZAK KH, RW MENDYK, and JJ WIENS. 2009. Can Parallel Diversification Occur in Sympatry? Repeated Patterns of Body-Size Evolution in Coexisting Clades of North American Salamanders. Evolution 63:1769-1784.}
		
		\item{WEISROCK DW, TJ PAPENFUSS, JR MACEY, SN LITVINCHUK, R POLYMENI, IH UGURTAS, E ZHAO, H JOWKAR, and A LARSON. 2006. A molecular assessment of phylogenetic relationships and lineage accumulation rates within the family Salamandridae (Amphibia, Caudata). Molecular Phylogenetics and Evolution 41:368-383.}
		
		\item{WIENS JJ and JT HOVERMAN. 2008. Digit reduction, body size, and paedomorphosis in salamanders. Evolution and Development 10:449-463.}
		
		\item{ZHANG P, Y-Q CHEN, H ZHOU, X-L WANG, TJ PAPENFUSS, DB WAKE and L-H QU. 2006. Phylogeny, evolution, and biogeography of Asiatic salamanders (Hynobiidae). Proceedings of the National Academy of Sciences USA 103:7360-7365.}
		
		\item{ZHANG P and DB WAKE. 2009. Higher-level salamander relationships and divergence dates inferred from complete mitochondrial genomes. Molecular Phylogenetics and Evolution. 53:492-508. }
}
}
