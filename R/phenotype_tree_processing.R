#'Converts a binary phenotype tree to a vector of foreground species names
#' @param tree a phylo object with binary edge lengths (foreground edges = 1, background edges = 0)
#' @return foregrounds a vector of foreground species names
#' @export
getForegroundsFromTree=function(tree){
  edge_length = tree$edge.length
  edge_nodes = tree$edge

  foreground_edges = edge_nodes[which(edge_length==1),]
  num_tips = length(tree$tip.label)

  tip_foregrounds = tree$tip.label[foreground_edges[which(foreground_edges[,2] <= num_tips),2]]
  ind_internal_foregrounds = which(foreground_edges[,2] > num_tips)
  if (length(ind_internal_foregrounds) > 0){
    internal_foregrounds = NULL
    for (i in 1:length(ind_internal_foregrounds)){
      foreground_branch = foreground_edges[ind_internal_foregrounds[i],]
      internal_foregrounds = c(internal_foregrounds, getInternalBranchName(foreground_branch[2], tree, foreground_edges, num_tips))
    }
    foregrounds = c(tip_foregrounds, internal_foregrounds)
  } else {
    foregrounds = tip_foregrounds
  }
  foregrounds
}

#'Produces phyloConverge-ready name of an ancestral branch
#' @param node internal node ID number, corresponding to the node ID names in the tree object
#' @param tree a phylo object with the tree topology of interest
#' @param allEdges a matrix containing all edges (e.g., the 'edge' component in the tree object, or a subset of it)
#' @param num_tips total number of tip species in the topology
#' @return branch_name the name of the input ancestral node
#' @export
getInternalBranchName=function(node,tree,allEdges,num_tips){
  daughter_nodes = getTipDaughterNodes(node,allEdges,num_tips)
  reprDaughters = tree$tip.label[c(min(daughter_nodes), max(daughter_nodes))]
  branch_name = paste(reprDaughters, collapse="-")
  branch_name
}

#' @keywords internal
getDaughterNodes=function(node, allEdges){
  daughters = allEdges[which(allEdges[,1] == node),2]
  daughters
}

#' @keywords internal
getTipDaughterNodes=function(node, allEdges, num_tips){
  if (node <= num_tips){
    daughters=NULL
  } else {
    alldaughters=allEdges[which(allEdges[,1]==node),2]
    ind_tip_daughters = which(alldaughters<=num_tips)
    if (length(ind_tip_daughters) == 2){
      daughters=alldaughters
    } else {
      ind_anc_daughters = which(alldaughters > num_tips)
      for (i in 1:length(ind_anc_daughters)){
        daughters = c(alldaughters[ind_tip_daughters],getTipDaughterNodes(alldaughters[ind_anc_daughters[i]], allEdges, num_tips))
      }
    }
  }
  daughters
}

#' @keywords internal
getTipDaughterNames=function(node,tree,allEdges,num_tips){
  daughter_nodes = getTipDaughterNodes(node,allEdges,num_tips)
  daughter_names = tree$tip.label[sort(daughter_nodes, decreasing=F)]
  daughter_names
}


#'Returns tree (phylo) object from evolution model
#' @param mod path to evolution model file
#' @return tree tree object in the evolution model
#' @export
getTreeFromEvolMod=function(mod){
  modObj = read.tm(mod)
  tree_nwk = modObj$tree
  tree = read.tree(text=tree_nwk)
  tree
}

#'Converts foreground species names into a binary phenotype tree with the given topology
#' @param foregrounds a character vector of foreground species
#' @param topology the intended tree topology represented as a phylo object
#' @param plotTree Boolean flag for plotting the resulting binary tree
#' @return treeout binary phenotype tree
#' @export
getTreeFromForegrounds=function(foregrounds, topology, plotTree=F){
  if (class(topology)=="phylo"){
    treeout = topology
  } else {
    stop('Wrong topology format.')
  }

  treeout$edge.length = rep(0, length(topology$edge.length))

  fgNodeIDs = unlist(lapply(foregrounds, findNodeID, tree=treeout))
  ind_fg_edges = which(treeout$edge[,2] %in% fgNodeIDs)
  treeout$edge.length[ind_fg_edges] = 1
  if (plotTree){
    plot(treeout)
  }
  treeout
}

#' @keywords internal
findNodeID=function(species_name, tree){
  if (grepl("-", species_name)){
    tip_daughters = strsplit(species_name, "-")[[1]]
    tip_daughters_ID = which(tree$tip.label %in% tip_daughters)
    nodeID = findCommonAncestorNode(tip_daughters_ID, tree)
  } else {
    nodeID = which(tree$tip.label == species_name)
  }
  nodeID
}

#' @keywords internal
findCommonAncestorNode=function(nodeIDs, tree){
  num_nodes = length(nodeIDs)

  allEdges = tree$edge

  destin_nodes = nodeIDs

  all_anc_nodes = NULL
  ancestor_found = FALSE
  while (!ancestor_found){
    ind_edges_iter = which(allEdges[,2] %in% destin_nodes)
    edges_iter = allEdges[ind_edges_iter,]

    source_nodes = edges_iter[,1]
    all_anc_nodes = c(all_anc_nodes, source_nodes)
    table_anc_nodes = table(all_anc_nodes)

    if (length(which(table_anc_nodes == length(nodeIDs)))){
      ancestor_found = TRUE
      ind_ca = which(table_anc_nodes == length(nodeIDs))
      ca = names(table_anc_nodes)[ind_ca]
    } else {
      destin_nodes = unique(source_nodes)
    }
  }

  as.numeric(ca)
}
