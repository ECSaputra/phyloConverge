require("rphast")

#'Subset foreground species based on existing species in alignment file
#' @param align alignment object read from read.msa
#' @param foregrounds character vector containing the foreground species
#' @return foregrounds.out character vector containing the edited foregrounds
#' @export
checkForegrounds=function(align, foregrounds){
  species.in.alignment = align$names
  tip.foregrounds = foregrounds[which(!grepl("-", foregrounds))]
  foregrounds.out = tip.foregrounds[which(tip.foregrounds %in% species.in.alignment)]

  internal.nodes = foregrounds[which(grepl("-", foregrounds))]
  if (length(internal.nodes) > 0){
    for (i in 1:length(internal.nodes)){
      internal.node = internal.nodes[i]
      daughter.tips = strsplit(internal.node,"-")[[1]]

      if (length(which(daughter.tips %in% species.in.alignment)) == 2){
        foregrounds.out = c(foregrounds.out, internal.node)
      }
    }
  }
  foregrounds.out
}


#'Identify which genetic elements in a table exist in a given multiple alignment file
#' @param elements_bed a BED format matrix/data frame containing a list of genetic elements with their chromosome and coordinates
#' @param maf an MSA object containing the sequence alignment
#' @return elements_in_maf a BED format matrix/data frame listing the elements that exist in maf
#' @export
getElementsInMaf=function(elements_bed, maf){
  coord_range = coord.range.msa(maf)

  ind_elements_in_maf = intersect(which(elements_bed[,2] >= coord_range[1]), which(elements_bed[,3] <= coord_range[2]))

  if (length(ind_elements_in_maf) > 0){
    elements_in_maf = elements_bed[ind_elements_in_maf,]
  } else {
    elements_in_maf = NULL
  }
  elements_in_maf
}

#'Convert a BED matrix/data.frame containing information on genetic elements into a features object
#' @param bed_file a BED format matrix/data frame containing a list of genetic elements with their chromosome, coordinates, and name (optional)
#' @param refseq the name of the sequence in the alignment that is used as a frame of reference
#' @param feature the feature type name (e.g., "exon", "intron", "CDS", etc. Default NULL)
#' @param strand a character string denoting the strand (e.g. "+", "-". Default NULL if strand is not relevant)
#' @param frame a 0, 1, or 2 specifying whether the feature is in frame
#' @param attribute a character string denoting the label of the feature
#' @return featureout a features object ready to use for phyloConverge
#' @export
convertBedToFeature=function(bed_file, refseq, feature=NULL, strand=NULL, frame=NULL, attribute=NULL){
  if (ncol(bed_file) < 3){
    stop("Incomplete bed file.")
  }

  if (is.null(attribute)){
    if (ncol(bed_file) > 3){
      attribute = bed_file[,4]
    } else {
      attribute = "."
    }
  }
  if (is.null(feature)){
    feature="."
  }
  featureout = feat(seqname= refseq, feature=".", start=bed_file[,2], end=bed_file[,3], attribute=attribute, strand=strand, frame=frame)
  featureout
}




#'Convert an entire alignment into a features object
#' @param aln msa object representing the alignment
#' @param refseq the name of the sequence in the alignment that is used as a frame of reference
#' @param feature the feature type name (e.g., "exon", "intron", "CDS", etc. Default NULL)
#' @param strand a character string denoting the strand (e.g. "+", "-". Default NULL if strand is not relevant)
#' @param frame a 0, 1, or 2 specifying whether the feature is in frame
#' @param attribute a character string denoting the label of the feature
#' @return featureout a features object ready to use for phyloConverge
#' @export
convertAlignmentToFeature=function(aln, refseq, feature=NULL, strand=NULL, frame=NULL, attribute=NULL){
  #aln = read.msa(alignment_file, refseq)
  aln_range = coord.range.msa(aln)
  if (is.null(feature)){
    feature='.'
  }
  if (is.null(attribute)){
    attribute='.'
  }
  featureout = feat(seqname=refseq, start=aln_range[1], end=aln_range[2], feature=feature, strand=strand, frame=frame, attribute=attribute)
  featureout
}
