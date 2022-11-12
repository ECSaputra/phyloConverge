
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phyloConverge

<!-- badges: start -->

<!-- badges: end -->

phyloConverge is a comparative genomics software for identifying
phenotypic associations of genetic elements.

phyloConverge combines maximum likelihood estimation and phylogeny-aware
permutation tests to detect significant convergent evolutionary rate
shifts underlying phenotypic convergence, while calibrating for biases
from phylogenetic dependence, sequence genomes, and correlations between
genetic elements. phyloConverge is able to analyze segments of
nucleotides at a wide range of sizes, from scoring genes or conserved
elements, to scanning entire multiple sequence alignments at nucleotide
resolution with computational tractability.

## Dependencies

phyloConverge requires the following to be installed:

  - *RERconverge*
    
      - Follow installation here:
        [RERconverge](https://github.com/nclark-lab/RERconverge)

  - RPHAST
    
      - Follow installation here:
        [RPHAST](https://github.com/CshlSiepelLab/RPHAST) or
        [here](https://confluence.cc.lehigh.edu/display/LKB/Install+RPHAST+on+Windows)

## Installation

phyloConverge can be installed using the devtools package, using the
following command:

``` r
library(devtools)
install_github("ECSaputra/phyloConverge")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(phyloConverge)
## basic example code
```
