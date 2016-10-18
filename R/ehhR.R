#####given a matrix, a function gives the ehh for each column from left to right.######
ehhLtoR=function(hapmatrix){
  require(plyr)
  hapmatrix=as.matrix(hapmatrix)
  den=count(hapmatrix[,1])[,2]*(count(hapmatrix[,1])[,2]-1)
  output=matrix(1,nrow=2,ncol=dim(hapmatrix)[2])
  if (dim(hapmatrix)[2]==1) return(c(1,1))
  ################
  if (dim(hapmatrix)[2]>1) {
    for (i in 2:dim(hapmatrix)[2]){
      cntmat=count(hapmatrix[,1:i])
      row1=which(cntmat[,1]==1)
      row2=which(cntmat[,1]==2)
      output[1,i]=sum(cntmat[row1,i+1]*(cntmat[row1,i+1]-1))/den[1]
      output[2,i]=sum(cntmat[row2,i+1]*(cntmat[row2,i+1]-1))/den[2]
    }
    return(output)
  }

}

#' A function to calculate EHH (Extended Haplotype Homozygosity) given a focal SNPs and a haplotype matrix.
#'
#' @param hapmatrix A haplotype matrix with ones (ancestral) and twos (derived).
#' @param site A specific site.
#' @return EHH (Extended Haplotype Homozygosity).
#' @examples
#' exmtx=matrix(sample(2,80,replace=T),ncol=10,nrow=8)
#' ehhR(exmtx,2)
#' ehhR(exmtx,4)
#'
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#'
#' @export
ehhR=function(hapmatrix,site){
  require(plyr)
  mend=dim(hapmatrix)[2]
  output=matrix(0,nrow=2,dim(hapmatrix)[2])
  output[,site]=c(1,1)
  output[,site:mend]=ehhLtoR(hapmatrix[,site:mend])
  output[,site:1]=ehhLtoR(hapmatrix[,site:1])
  rownames(output)=c("Anc. Allele","Der. Allele")
  return(output)
}
