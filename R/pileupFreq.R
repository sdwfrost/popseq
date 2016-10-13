#' Converts a pileup to a dataframe of nucleotide frequencies
#'
#' \code{pileupFreq} converts the output of \code{pileup} to a dataframe of nucleotide frequencies.
#'
#' @param p a pileup
#' 
#' @return Nucleotide frequencies (as a data frame)
#'
#' @author Simon Frost (\email{sdwfrost@@gmail.com})
#'
#' @export
pileupFreq <- function(p) {
  nucleotides <- levels(p$nucleotide)
  res <- split(p, p$seqnames)
  res <- lapply(res, function (x) {split(x, x$pos)})
  res <- lapply(res, function (positionsplit) {
    nuctab <- lapply(positionsplit, function(each) {
      chr = as.character(unique(each$seqnames))
      pos = as.character(unique(each$pos))
      tablecounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
      c(chr,pos, tablecounts)
    })
    nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
    rownames(nuctab) <- NULL
    nuctab
  })
  res <- data.frame(do.call("rbind", res),stringsAsFactors=F)
  rownames(res) <- NULL
  colnames(res) <- c("seqnames","start",levels(p$nucleotide))
  res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
  res
}
