##############hamming distance matrix############
hammatrix=function(freqMatrix){
  require(stringdist)
  hammatrix=matrix(0,nrow=dim(freqMatrix)[1],ncol=dim(freqMatrix)[1])
  rownames(hammatrix)=rownames(freqMatrix)
  colnames(hammatrix)=rownames(freqMatrix)
  for (i in 1:dim(freqMatrix)[1]){
    for (j in 1:dim(freqMatrix)[1]){
      hammatrix[i,j]=stringdist(rownames(freqMatrix)[i],rownames(freqMatrix)[j],method="hamming")
    }
  }
  return(hammatrix)
}

####the mutation rate matrix############
mutationcoe=function(freqMatrix,rate){
  L=length(rownames(freqMatrix))/dim(freqMatrix)[1]
  #######the hamming distance between two haplotype sequences######
  #rate=1/3*10^{-5}
  require(stringdist)
  mu=matrix(0,nrow=dim(freqMatrix)[1],ncol=dim(freqMatrix)[1])
  rownames(mu)=rownames(freqMatrix)
  colnames(mu)=rownames(freqMatrix)
  for (i in 1:dim(freqMatrix)[1]){
    for (j in 1:dim(freqMatrix)[1]){
      mu[i,j]=rate^(stringdist(rownames(mu)[i], rownames(mu)[j],method="hamming"))*(1-3*rate)^(L-stringdist(rownames(mu)[i], rownames(mu)[j],method="hamming"))
    }
  }
  return(mu)
}


##########selection effect on the largest minor allel########
fnselectsingle=function(frequencyvec,mutationMatrix,selecoe){
  ###frequencyvec is a frequency vector of 4 elements; selecoe is a selection coefficient vector of 4 elements, the non-zero one is the selection coefficient;###
  ###when this is used for a neutral model the selecoe==c(0,0,0,0); mutationMatrix is 4*4 matrix of mutation coefficient.##############
  #newfreq=0
  #newfreq[1]=frequencyvec[1]*exp(log(selecoe)*2)/(frequencyvec[1]*exp(log(selecoe)*2)+1-frequencyvec[1])
  #newfreq[2]=1-newfreq[1]

  newfreq=mutationMatrix%*%mutationMatrix%*%frequencyvec
  newfreq=newfreq*selecoe
  newfreq=newfreq/sum(newfreq)
  return(newfreq)
}

################tilde frequency function##########
fnqtilde=function(frequencyvec,hammingmatrix,e){
  newvec=0
  for (i in 1:length(frequencyvec)){
    newvec[i]=frequencyvec[i]+sum((frequencyvec[which(hammingmatrix[i,]==1)]-frequencyvec[i])*e)
  }
  return(newvec)
}



#' A function to calculate loglikelihood of a single locus model.
#' @param parr Parameter vector.
#' @param freqmatrix Longitudinal frequency matrix of 4 nucleotides (A, C, G, T).
#' @param Pos The position of A, C, G, T under selection.
#' @return loglikelihood of a single locus model.
#' @examples
#' freq=cbind(D2=c(0,0,57,0),D3=c(0,0,44,0),D4=c(3,0,45,0))
#' rownames(freq)=c("a","c","g","t")
#' flogliksingle(c(0,1,1,1,1,1),freq,1/3*10^(-5),1)
#' library(minqa)
#' bobyqa(c(.1,.1,.1,.1,.1,.1),flogliksingle,upper=Inf,lower=0,freqMatrix=freq,rate=1/3*10^(-5),pos=1)
#'
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#'
#' @export
#'
flogliksingle=function(parr,freqMatrix,rate,pos){
  ###freqmatrix is the frequency matrix;###
  ###rate is the mutation rate;############
  ###pos is the postion under selection;###
  ###par is the parameter vector to optimise###
  freqlth=dim(freqMatrix)[1]
  ham=hammatrix(freqMatrix)
  mu=mutationcoe(freqMatrix,rate)
  e=parr[1]
  inifreq=c(parr[2],parr[3],parr[4],parr[5])
  coevec=rep(1,4)
  if (pos != 0) {coevec[pos]=parr[6]}
  ####the no. of time units###
  T=dim(freqMatrix)[2]
  #####frequncy vector q givent the initial frequcency####/sum(inifreq)#
  freqmatrix=matrix(0,nrow=length(inifreq),ncol=T)
  ###from t0 to t1####
  q1=inifreq/sum(inifreq)
  freqmatrix[,1]=q1
  for (i in 2:T){
    q=fnselectsingle(freqmatrix[,i-1],mu,coevec)
    q=fnselectsingle(q,mu,coevec)
    freqmatrix[,i]=q
  }




  ####then the tilde freqmatrix#####
  tildfreqmatrix=matrix(0,nrow=length(inifreq),ncol=T)
  for (i in 1:T){
    tildfreqmatrix[,i]=fnqtilde(freqmatrix[,i],ham,e)
  }
  ######loglikelihood##########
  loglik=0
  for (i in 1:T){
    loglik[i]=dmultinom(freqMatrix[,i],size=sum(freqMatrix[,i]),prob=tildfreqmatrix[,i],log=T)
  }
  return(-sum(loglik))
}
