####a function to calculate the Linkage disequilibrium given two columns of 1s and 2s####
ld2=function(twocol){
  require(plyr)
  freq1=count(twocol[,1])
  freq2=count(twocol[,2])
  freq12=count(c(twocol[,1],twocol[,2]))
  p1=freq1[2,2]/sum(freq1[,2])
  p2=freq2[2,2]/sum(freq2[,2])
  p12=freq12[2,2]/sum(freq12[,2])
  num=p12-p1*p2
  den=p1*p2*(1-p1-p2+p1*p2)
  return(num^2/den)
}

####a function to calculate the Linkage disequilibrium of a sequences*SNPs matrix####
ldall=function(hapmatrix){
  require(plyr)
  hapmatrix=as.matrix(hapmatrix)
  output=matrix(1,nrow=dim(hapmatrix)[2],ncol=dim(hapmatrix)[2])
  for (i in 1:(dim(hapmatrix)[2]-1)){
    for (j in (i+1):dim(hapmatrix)[2]){
      output[i,j]=ld2(cbind(hapmatrix[,i],hapmatrix[,j]))
      output[j,i]=output[i,j]
    }
  }
  return(output)
}



#' A function to calculate Kelly's Zns statistic and the omega_max statistic given a haplotype matrix.
#'
#' @param hapmatrix A haplotype matrix with ones (ancestral) and twos (derived).
#' @return Kelly's Zns statistic, the omega_max statistic and the site produces omega_max.
#' @examples
#' exmtx=matrix(sample(2,80,replace=T),ncol=10,nrow=8)
#' ldR(exmtx)
#'
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#'
#' @export
#'
ldR=function(hapmatrix){
  require(plyr)
  rmat=ldall(hapmatrix)
  S=dim(hapmatrix)[2]
  zns=(sum(rmat)-S)/S/(S-1)
  omega=0
  for (i in 1:(S-1)){
    omega[i]=(sum(rmat[1:i,1:i])+sum(rmat[(i+1):S,(i+1):S])-S)/2/(choose(i,2)+choose(S-i,2))
    omega[i]=omega[i]*i*(S-i)/(sum(rmat[1:i,(i+1):S]))
  }

  return(list(Zns=zns,omega_max=c(max(omega),which(omega==max(omega)))))
}
