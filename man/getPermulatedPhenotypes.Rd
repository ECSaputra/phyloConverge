% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permulation_functions.R
\name{getPermulatedPhenotypes}
\alias{getPermulatedPhenotypes}
\title{Produce permulated phenotype trees from the observed phenotype and phylogeny}
\usage{
getPermulatedPhenotypes(
  foregrounds,
  neutraltree,
  num_perms,
  root_species,
  output_mod = "names"
)
}
\arguments{
\item{foregrounds}{a vector of foreground species names}

\item{neutraltree}{a phylo object representing neutral evolution}

\item{num_perms}{number of permulations}

\item{root_species}{the species to root the tree on}

\item{output_mod}{flag for the format of the output ("names" outputs permulated foreground names, "tree" outputs permulated phenotype trees)}
}
\value{
a list object containing permulated phenotypes
}
\description{
Produce permulated phenotype trees from the observed phenotype and phylogeny
}
