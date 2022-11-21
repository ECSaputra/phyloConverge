require('RERconverge')

#'Produce a permulated phenotype tree from the observed phenotype and phylogeny
#' @param foregrounds a vector of foreground species names
#' @param neutraltree a phylo object representing neutral evolution
#' @param root_species the species to root the tree on
#' @return permulated_tree a permulated phenotype tree with binary edge lengths (foreground edges = 1, background edges = 0)
#' @export
getPermulatedPhenotypeTree=function(foregrounds, neutraltree, root_species){
  foreground_tree = getTreeFromForegrounds(foregrounds, neutraltree, plotTree=F)
  binary_permulation_inputs = getBinaryPermulationInputsFromTree(foreground_tree)
  sisters_list = binary_permulation_inputs$sisters_list
  fg_vec = binary_permulation_inputs$fg_vec

  pathvec = getPathvecfromForegroundTree(foreground_tree)

  masterTree = list()
  masterTree[[1]] = neutraltree
  names(masterTree) = c("masterTree")

  permulated_tree = simBinPhenoCC(masterTree, neutraltree, root_species, fg_vec, sisters_list=sisters_list, pathvec, plotTreeBool = F)

  permulated_tree
}

#'Produce permulated phenotype trees from the observed phenotype and phylogeny
#' @param foregrounds a vector of foreground species names
#' @param neutraltree a phylo object representing neutral evolution
#' @param num_perms number of permulations
#' @param root_species the species to root the tree on
#' @param output_mod flag for the format of the output ("names" outputs permulated foreground names, "tree" outputs permulated phenotype trees)
#' @return a list object containing permulated phenotypes
#' @export
getPermulatedPhenotypes=function(foregrounds, neutraltree, num_perms, root_species, output_mod="names"){
  fg_reps = lapply(1:num_perms, return_object, x=foregrounds)
  permulated_trees = lapply(fg_reps, getPermulatedPhenotypeTree, neutraltree=neutraltree, root_species=root_species)

  if (output_mod == "trees"){
    return(permulated_trees)
  } else if (output_mod == "names"){
    permulated_foregrounds = lapply(permulated_trees, getForegroundsFromTree)
    return(permulated_foregrounds)
  }
}

#' @keywords internal
getPathvecfromForegroundTree=function(foreground_tree){
  fg_edge_idx = which(foreground_tree$edge.length==1)
  fg_node_idx = foreground_tree$edge[fg_edge_idx,2]
  idx_tip_fg = which(fg_node_idx <= length(foreground_tree$tip.label))
  foreground_species = foreground_tree$tip.label[fg_node_idx[idx_tip_fg]]

  pathvec = rep(0, length(foreground_tree$tip.label))
  names(pathvec) = foreground_tree$tip.label
  pathvec[foreground_species] = 1
  pathvec
}

#' @keywords internal
return_object=function(x_idx, x){
  return(x)
}

