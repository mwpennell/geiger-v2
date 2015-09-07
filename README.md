# Geiger

[![Build Status](https://travis-ci.org/mwpennell/geiger-v2.png?branch=master)](https://travis-ci.org/mwpennell/geiger-v2)

This is the development version of the [geiger package](http://bioinformatics.oxfordjournals.org/content/30/15/2216) for manipulating phylogenetic comparative data
and fitting macroevolutionary models. The package can be downloaded from CRAN
```
install.packages("geiger")
```
or installed directly from github with the [devtools](https://github.com/hadley/devtools) package
```
install.packages("devtools")
install_github("mwpennell/geiger-v2")
```

## Major features

geiger is a (growing) collection of methods developed over the years by many researchers.
Here is a a non-comprehensive list of methods:

* Fit continuous models of evolution (BM, OU, EB, Pagel models, etc.)
* Fit discrete models of evolution (Mk and variants)
* [Identify shifts in the rate of continuous trait evolution](http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2011.01401.x/abstract)
* [Fit continuous trait models to unresolved data using ABC](http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2011.01474.x/abstract)
* [Use fossil information to improve macroevolutionary inference](http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2012.01723.x/abstract)
* [Identify shifts in the rate of diversification](http://www.pnas.org/content/106/32/13410.short)
* [Posterior predictive model assessment](http://sysbio.oxfordjournals.org/content/63/3/293)
* [Time-scaling large phylogenies with 'congruification'](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12051/full)

## Citing geiger

If you use geiger, please cite:

Pennell, M.W., J.M. Eastman, G.J. Slater, J.W. Brown, J.C. Uyeda, R.G. FitzJohn, M.E. Alfaro, and L.J. Harmon. 2014. geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. Bioinformatics 15:2216-2218.

in addition to the original papers describing the methods.

## Feedback

We are always looking to improve geiger. If you have comments/questions/ideas, we encourage you to get in contact by posting an issue or making a pull request.
