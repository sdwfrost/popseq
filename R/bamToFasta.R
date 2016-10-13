#' Writes a BAM file to a file in FASTA format
#'
#' \code{bamToFasta} converts a BAM file to a FASTA file.
#'
#' @param bf a \code{BamFile}
#' @param fn the filename for the FASTA file
#' @param region an optional region from the BAM file (as a \code{GRanges} object)
#'
#' @author Simon Frost (\email{sdwfrost@@gmail.com})
#'
#' @export
bamToFasta <- function(bf,fn,region=NULL){
  if(!is.null(region)){
    s <- GenomicAlignments::stackStringsFromBam(bf, param=region,D.letter="-",N.letter="-", Lpadding.letter="-",Rpadding.letter="-",use.names=TRUE)
  }
  else{
    s <- GenomicAlignments::stackStringsFromBam(bf, D.letter="-",N.letter="-", Lpadding.letter="-",Rpadding.letter="-",use.names=TRUE)
  }
  Biostrings::writeXStringSet(s,file=fn,format="fasta")
}
