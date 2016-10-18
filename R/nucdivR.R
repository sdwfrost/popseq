
#' A function to calculate the average hamming distance within each two populations, between two populations and Fst statistic.
#' @param hapmatrix1 A haplotype matrix with ones (ancestral) and twos (derived).
#' @param hapmatrix2 The other haplotype matrix with ones (ancestral) and twos (derived).
#' @return The average hamming distance within each two populations, between two populations and Fst statistic.
#' @examples
#' exmtx=matrix(sample(2,80,replace=T),ncol=10,nrow=8)
#' exmtx1=matrix(sample(2,80,replace=T),ncol=10,nrow=8)
#' nucdivR(exmtx,3-exmtx)
#' nucdivR(exmtx,exmtx1)
#'
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#'
#' @export
nucdivR=function(hapmatrix1, hapmatrix2){
  require(e1071)
  hapmatrix1=as.matrix(hapmatrix1)
  hapmatrix2=as.matrix(hapmatrix2)
  S1=dim(hapmatrix1)[2]
  S2=dim(hapmatrix2)[2]
  pi1=sum(hamming.distance(t(hapmatrix1)))/(S1-1)^2
  pi2=sum(hamming.distance(t(hapmatrix2)))/(S2-1)^2
  hap12=cbind(hapmatrix1,hapmatrix2)
  pi12=(sum(hamming.distance(t(hap12)))-sum(hamming.distance(t(hapmatrix1)))-sum(hamming.distance(t(hapmatrix2))))/S1/S2/2
  return(list(pi_x=pi1,pi_y=pi2,d_xy=pi12,Fst_XY=1-(pi1+pi2)/2/pi12))
}
