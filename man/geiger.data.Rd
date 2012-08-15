\name{geiger.data}
\alias{amphibia}
\alias{chelonia}
\alias{geospiza}
\alias{primates}
\alias{ncbi}
\alias{caudata}

\title{example datasets}
\description{Provides several datasets and trees}
\usage{data(amphibia); data(caudata); data(chelonia); data(geospiza); data(primates); data(ncbi)} 

\details{
The objects \code{caudata}, \code{chelonia}, \code{geospiza}, and \code{primates} each have at least two items, a \code{phy} object 
and a \code{dat} object. The \code{phy} object is a phylogenetic tree of class 'phylo' (see \code{\link[ape]{read.tree}}).  The \code{dat} object 
is a named vector of (natural log-transformed) body sizes for each group. The salamander data (\code{caudata}) also includes 
the systematics for all recognized taxa (\code{tax}) as well as a time-calibrated family-level phylogeny (\code{fam})from Zhang and Wake (2009).

The \code{amphibia} object is a set of three trees, the last of which is a time-scaled estimate of the Pyron and Wiens (2011) tree (see \code{\link[geiger]{congruify.phylo}}).

The \code{ncbi} object is a copy of the NCBI taxonomy resource (see \url{http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html}). 
}

\references{

Data are from the following sources:

\code{AMPHIBIANS (amphibia)}

\itemize{
		\item{ROELANTS K, DJ GOWER, M WILKINSON, SP LOADER, SD BIJU, K GUILLAUME, L MORIAU, and F BOSSUYT. 2007. Global patterns of diversification in the history of modern amphibians. PNAS 104:887-892.}
		
		\item{PYRON RA and JJ WIENS. 2011. A large-scale phylogeny of Amphibia including over 2,800 species, and a revised classification of extant frogs, salamanders, and caecilians. Molecular Phylogenetics and Evolution 61:543–583.}
}

\code{SALAMANDERS (caudata)}

\itemize{

		\item{ADAMS DC, CM BERNS, KH KOZAK, and JJ WIENS. 2009. Are rates of species diversification correlated with rates of morphological evolution? Proceedings of the Royal Society of London B 276:2729-2738. }
		
		\item{BONETT RM, PT CHIPPINDALE, PE MOLER, RW VAN DEVENDER, and DB WAKE. 2009. Evolution of gigantism in amphiumid salamanders. PLoSONE 4(5):e5615.}
				
		\item{KOZAK KH, RW MENDYK, and JJ WIENS. 2009. Can Parallel Diversification Occur in Sympatry? Repeated Patterns of Body-Size Evolution in Coexisting Clades of North American Salamanders. Evolution 63:1769-1784.}
		
		\item{WEISROCK DW, TJ PAPENFUSS, JR MACEY, SN LITVINCHUK, R POLYMENI, IH UGURTAS, E ZHAO, H JOWKAR, and A LARSON. 2006. A molecular assessment of phylogenetic relationships and lineage accumulation rates within the family Salamandridae (Amphibia, Caudata). Molecular Phylogenetics and Evolution 41:368-383.}
		
		\item{WIENS JJ and JT HOVERMAN. 2008. Digit reduction, body size, and paedomorphosis in salamanders. Evolution and Development 10:449-463.}
		
		\item{ZHANG P, Y-Q CHEN, H ZHOU, X-L WANG, TJ PAPENFUSS, DB WAKE and L-H QU. 2006. Phylogeny, evolution, and biogeography of Asiatic salamanders (Hynobiidae). Proceedings of the National Academy of Sciences USA 103:7360-7365.}
		
		\item{ZHANG P and DB WAKE. 2009. Higher-level salamander relationships and divergence dates inferred from complete mitochondrial genomes. Molecular Phylogenetics and Evolution. 53:492-508. }
}

\code{PRIMATES (primates)}

\itemize{
		\item{REDDING DW, C DEWOLFF, and AO MOOERS. 2010. Evolutionary distinctiveness, threat status and ecological oddity in primates. Conservation Biology 24:1052-1058.}
		
		\item{VOS RA and AO MOOERS. 2006. A new dated supertree of the Primates. Chapter 5. In: VOS RA (Ed.) Inferring large phylogenies: the big tree problem. [Ph.D. thesis]. Burnaby BC, Canada: Simon Fraser University.}
}
	

	
\code{TURTLES (chelonia)}

\itemize{
		\item{JAFFE, G SLATER, and M ALFARO. 2011. Ecological habitat and body size evolution in turtles. Biology Letters: in press. }
}	


\code{DARWINS FINCHES (geospiza)}

\itemize{

		\item{SCHLUTER D}

}


\code{NCBI}

\itemize{
		\item{SAYERS EW, T BARRETT, DA BENSON, SH BRYANT, K CANESE, V CHETVERNIN, DM CHURCH, M DICUCCIO, R EDGAR, S FEDERHEN, M FEOLO, LY GEER, M HLEMBERG, Y KAPUSTIN, D LANDSMAN, DJ LIPMAN, TL MADDEN, DR MAGLOTT, V MILLER, I MIZRACHI, J OSTELL, KD PRUITT, GD SCHULER, E SEQUEIRA, ST SHERRY, M SHUMWAY, K SIROTKIN, A SOUVOROV, G STARCHENKO, TA TATUSOVA, L WAGNER, E YASCHENKO, and J YE. 2009. Database resources of the National Center for Biotechnology Information. Nucleic Acids Research 37:D5-15.}
   
		\item{BENSON DA, I KARSCH-MIZRACHI, DJ LIPMAN, J OSTELL, and EW SAYERS. 2009. GenBank. Nucleic Acids Research 37:D26-31.}
}



}
