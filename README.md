# Geiger

[![Build Status](https://travis-ci.org/mwpennell/geiger-v2.svg?branch=master)](https://travis-ci.org/mwpennell/geiger-v2)

This is the development version of the [geiger package](https://academic.oup.com/bioinformatics/article/30/15/2216/2390619) for manipulating phylogenetic comparative data
and fitting macroevolutionary models. The package can be downloaded from CRAN
```
install.packages("geiger")
```
or installed directly from github with the [devtools](https://github.com/r-lib/devtools) package
```
install.packages("devtools")
library(devtools)
install_github("mwpennell/geiger-v2")
```

## Major features

geiger is a collection of methods developed over the years by many researchers.
Here is a a non-comprehensive list of methods:

* Fit continuous models of evolution (BM, OU, EB, Pagel models, etc.)
* Fit discrete models of evolution (Mk and variants)
* Identify shifts in the rate of continuous trait evolution
* Fit continuous trait models to unresolved data using ABC
* Use fossil information to improve macroevolutionary inference
* Identify shifts in the rate of diversification
* Posterior predictive model assessment
* Time-scaling large phylogenies with 'congruification'

## Citing geiger

If you use geiger, please cite:

Pennell, M.W., J.M. Eastman, G.J. Slater, J.W. Brown, J.C. Uyeda, R.G. FitzJohn, M.E. Alfaro, and L.J. Harmon. 2014. geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. Bioinformatics 30:2216-2218.

in addition to the original papers describing the methods.

## Acknowledgements

We thank the CRAN team for help cleaning up our package errors.

## Feedback

We are always looking to improve geiger. If you have comments/questions/ideas, we encourage you to get in contact by posting an issue or making a pull request.
