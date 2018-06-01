\name{r8s.phylo}
\alias{r8s.phylo}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{
call r8s from geiger}
\description{
call r8s, including a calibration file}
\usage{
r8s.phylo(phy, calibrations=NULL, base="r8srun", ez.run="none", rm=TRUE,
  blformat=c(lengths="persite", nsites=10000, ultrametric="no", round="yes"),
  divtime=c(method="NPRS", algorithm="POWELL"),
  cv=c(cvStart=0, cvInc=0.5, cvNum=8), do.cv=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogram to turn into a chronogram}
  \item{calibrations}{a set of calibrations}
  \item{base}{file name}
  \item{ez.run}{if set to "PL" or "NPRS", does the analysis with settings appropriate for each}
  \item{rm}{remove the output files}
  \item{blformat}{blformat options for r8s. Pay special attention to nsites.}
  \item{divtime}{divtime options for r8s}
  \item{cv}{cv options for r8s}
  \item{do.cv}{Boolean for whether to do cross validation or not}
}
\details{
	This function uses r8s and a calibration to make your tree ultrametric}

\references{
SANDERSON r8s ____NEED TO ADD_____
}

\author{JM Eastman & B O'Meara}

\examples{

\dontrun{
phy <- read.tree(text=paste0("(Marchantia:0.033817,",
  "(Lycopodium:0.040281,((Equisetum:0.048533",
  "Osmunda:0.033640,Asplenium:0.036526):0.000425):",
  "0.011806,((((Cycas:0.009460,Zamia:0.018847):",
  "0.005021,Ginkgo:0.014702):1.687e-86,((Pinus:",
  "0.021500,(Podocarpac:0.015649,Taxus:0.021081):",
  "0.006473):0.002448,(Ephedra:0.029965,(Welwitsch",
  ":0.011298,Gnetum:0.014165):0.006883):0.016663)",
  ":0.006309):0.010855,((Nymphaea:0.016835,(((((Saururus:",
  "0.019902,Chloranth:0.020151):1.687e-86,",
  "((Araceae:0.020003,(Palmae:0.006005,Oryza:0.031555):",
  "0.002933):0.007654,Acorus:0.038488):0.007844)",
  ":1.777e-83,(Calycanth:0.013524,Lauraceae:0.035902):",
  "0.004656):1.687e-86,((Magnolia:0.015119,Drimys:",
  "0.010172):0.005117,(Ranunculus:0.029027,((Nelumbo:",
  "0.006180,Platanus:0.002347):0.003958,(Buxaceae:",
  "0.013294,((Pisum:0.035675,(Fagus:0.009848,Carya:",
  "0.008236):0.001459):0.001994,(Ericaceae:0.019136,",
  "Solanaceae:0.041396):0.002619):1.687e-86):0.004803)",
  ":1.687e-86):0.006457):0.002918):0.007348,",
  "Austrobail:0.019265):1.687e-86):1.687e-86,Amborella:",
  "0.019263):0.003527):0.021625):0.012469):",
  "0.019372);"))

	calibrations <- data.frame(MRCA="LP", MaxAge=450, MinAge=450,
		taxonA="marchantia", taxonB="pisum", stringsAsFactors=FALSE)

	phy.nprs <- r8s.phylo(phy=phy, calibrations=calibrations, base="nprs_file", ez.run="NPRS")
	phy.pl <- r8s.phylo(phy=phy, calibrations=calibrations, base="pl_file", ez.run="PL")
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphs}
\keyword{manip}
