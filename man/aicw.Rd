\name{aicw}
\alias{aicw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
determining Akaike weights
}
\description{
This function computes Akaike Weights and ranks model based on their support from a vector of AIC scores. 
}
\usage{
aicw(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a named vector of AIC scores
}
}
\value{
%%  ~Describe the value returned
An nx3 matrix, where n is the number of models being compared. The first column contains the AIC scores, the second contains the deltaAIC score and the third the Akaike Weight. Models are ranked in descending order according to weight. Rownames are the model names taken from the input vector.
}
\author{
Graham Slater
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
AIC.scores <- c(3,7,-5, 6)
names(AIC.scores) <- c("model1", "model2", "model3", "model4")
aicw(AIC.scores)
}
