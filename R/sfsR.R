#' A function to calculate Tajima's D and Fay and WU's H statsitc given a haplotype matrix.
#' @param hapmatrix A haplotype matrix with ones (ancestral) and twos (derived).
#' @return Tajima's D and Fay and WU's H statsitc.
#' @examples
#' exmtx=matrix(sample(2,80,replace=T),ncol=10,nrow=8)
#' sfsR(exmtx)
#' sfsR(3-exmtx)
#'
#'
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#'
#' @export
#'
sfsR=function(hapmatrix){
  n=dim(hapmatrix)[1]
  S=dim(hapmatrix)[2]
  I=1:(n-1)
  a1=sum(1/I)
  a2=sum(1/I^2)
  b1=(n+1)/3/(n-1)
  b2=2*(n^2+n+3)/9/n/(n-1)
  c1=b1-1/a1
  c2=b2-(n+2)/a1/n+a2/a1^2
  e1=c1/a1
  e2=c2/(a1^2+a2)
  Dnum=sum(hamming.distance(hapmatrix))/n/(n-1)-S/a1
  Dden=sqrt(e1*S+e2*S*(S-1))
  i=0
  for (j in 1:S){
    i[j]=length(which(hapmatrix[,j]!=hapmatrix[1,j]))
  }

  ii=as.numeric(names(table(i)))
  zeta=as.numeric(table(i))
  g1=(n-2)/6/(n-1)
  g2=(18*n^2*(3*n+2)*sum(1/(1:n)^2)-88*n^3-9*n^2+13*n-6)/9/n/(n-1)^2
  theta=S/a1
  theta2=S*(S-1)/(a1^2+b1)
  Hnum=sum(hamming.distance(hapmatrix))/n/(n-1)-sum(ii*zeta)/(n-1)
  Hden=sqrt(g1*theta+g2*theta2)

  return(list(TajimaD=Dnum/Dden,FayandWuH=Hnum/Hden))

}
