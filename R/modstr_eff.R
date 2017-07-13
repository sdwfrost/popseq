logPDM_eff=function(z,pstar=.1){
  ####z is a frequency vector of length 4, pstar is a mutation rate###
  p=0
  S=sum(z)
  comn1=pbeta(pstar,S-z[4]+1,z[4]+1,log.p=T)+lbeta(S-z[4]+1,z[4]+1)-log(pstar)
  #comn2=log(1-pbeta(pstar,S-z[4]+1,z[4]+1))+lbeta(S-z[4]+1,z[4]+1)-log(1-pstar)
  comn2=pbeta(pstar,S-z[4]+1,z[4]+1,log.p=T,lower.tail = F)+lbeta(S-z[4]+1,z[4]+1)-log(1-pstar)
  p[1]=(S-z[4])*log(1/3)+comn1
  p[2]=(z[2]+z[3])*log(1/2)+lfactorial(z[1])+lfactorial(z[2]+z[3])-lfactorial(S-z[4]+1)+comn1
  p[3]=(z[1]+z[3])*log(1/2)+lfactorial(z[2])+lfactorial(z[1]+z[3])-lfactorial(S-z[4]+1)+comn1
  p[4]=(z[1]+z[2])*log(1/2)+lfactorial(z[3])+lfactorial(z[1]+z[2])-lfactorial(S-z[4]+1)+comn1
  p[5]=log(2)+lfactorial(z[1])+lfactorial(z[2])+lfactorial(z[3])-lfactorial(S-z[4]+2)+comn1
  p[6]=(S-z[4])*log(1/3)+comn2
  p[7]=(z[2]+z[3])*log(1/2)+lfactorial(z[1])+lfactorial(z[2]+z[3])-lfactorial(S-z[4]+1)+comn2
  p[8]=(z[1]+z[3])*log(1/2)+lfactorial(z[2])+lfactorial(z[1]+z[3])-lfactorial(S-z[4]+1)+comn2
  p[9]=(z[1]+z[2])*log(1/2)+lfactorial(z[3])+lfactorial(z[1]+z[2])-lfactorial(S-z[4]+1)+comn2
  p[10]=log(2)+lfactorial(z[1])+lfactorial(z[2])+lfactorial(z[3])-lfactorial(S-z[4]+2)+comn2
  return(p)
}



#####an efficient ppa generation code given frequency matrix, pstar, priorPA and model structure matrix####
ppa.gen_eff=function(freqmatrix,pstar=.1,priorPA=.001,models,delta){
  ####delta is the indicator of models of interests
  nsamp=dim(models)[2]/2
  xx=models
  xx.index=unique(xx[,-1:-nsamp])
  ###seperate xx in to a list of matrix by the unque model indicator
  mat.list=list()
  for (i in 1:dim(xx.index)[1]){
    mat.list[[i]]=xx[apply(xx,1,function(x) all(x[-1:-nsamp]==xx.index[i,])),][,1:nsamp]
    mat.list[[i]]=rbind(xx.index[i,],mat.list[[i]])
  }
  #####do the same procedures for all separated matrix in the list #####
  logPDM.list=lapply(mat.list,function(x,y,pstar){
    df.mod =data.frame(Index=x[1,],mod=t(x[-1,]))
    df.freq=data.frame(Index=x[1,],freq=t(y))
    df.freq.merge=ddply(df.freq,"Index",numcolwise(sum))
    df.mod.merge=ddply(df.mod,"Index",unique)
    freqmat.merge=df.freq.merge[,-1]
    mod.merge=df.mod.merge[,-1]
    logfull=apply(as.matrix(freqmat.merge),1,logPDM_eff,pstar=pstar)
    lognmod=rbind(logfull,t(mod.merge)+1)
    result=apply(lognmod,2,function(x) x[1:10][x[-1:-10]])
    rowSums(result)


  }, y=freqmatrix,pstar=pstar)


  ########vector for p(m_k)####
  pmk=delta*priorPA/sum(delta)+(1-delta)*(1-priorPA)/sum(1-delta)
  result.PDM=do.call("c",logPDM.list)
  logPPA=result.PDM+log(pmk)
  ###log-sum-exp approximation
  lse.PPA=max(logPPA)+log(sum(exp(logPPA-max(logPPA))))
  PPA=exp(logPPA-lse.PPA)
  ###PPA=PPA/sum(PPA)

  ###############output#####################
  return(list(Model.Top5=models[order(-PPA),][1:5,],
              PPA.Top5=PPA[order(-PPA)][1:5],
              PPA_SI=sum(PPA*delta) ))


  }
