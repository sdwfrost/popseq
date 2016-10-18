######given a matrix, a function gives the ehhs for each column from left to right.####
ehhsLtoR=function(hapmatrix){
  require(plyr)
  hapmatrix=as.matrix(hapmatrix)
  den=sum(count(hapmatrix[,1])[,2])-sum((count(hapmatrix[,1])[,2])^2)
  output=1
  if (dim(hapmatrix)[2]==1) return(1)
  ################
  if (dim(hapmatrix)[2]>1) {
    for (i in 2:dim(hapmatrix)[2]){
      cntmat=count(hapmatrix[,1:i])
      nume=sum(count(hapmatrix[,1])[,2])-sum(cntmat[,i+1]^2)
      output[i]=nume/den
    }
    return(output)
  }

}

#' A function to calculate  EHHs (site Extended Haplotype Homozygosity) given a focal site and a haplotype matrix.
#'
#' @param hapmatrix A haplotype matrix with ones (ancestral) and twos (derived).
#' @param site A specific site.
#' @return EHH (Extended Haplotype Homozygosity).
#' @examples
#' exmtx=matrix(sample(2,80,replace=T),ncol=10,nrow=8)
#' ehhsR(exmtx,2)
#' ehhsR(exmtx,4)
#'
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#'
#' @export
ehhsR=function(hapmatrix,site){
  require(plyr)
  mend=dim(hapmatrix)[2]
  output=rep(0,mend)
  output[site:mend]=ehhsLtoR(hapmatrix[,site:mend])
  output[site:1]=ehhsLtoR(hapmatrix[,site:1])
  return(output)
}
