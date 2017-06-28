#' @useDynLib popseq genmodels
modstr<-function(nsamp1) .Call(genmodels,nsamp1)

###function to generate all possible models given a sample size####
modmat=function(sample.size){
models_num <- modstr(sample.size)
models_num <- matrix(models_num, ncol = 2 * sample.size)
#unique(models_num)
###indicator vector given the less stringent criteria where 1 indicate a model of interest
delta=rep(0,dim(models_num)[1])
indi=which(models_num[,2]>4 & models_num[,(2+sample.size)]>0)

if (sample.size>2){
  for (i in 3:sample.size){
    d=which(models_num[,i]>4 & models_num[,(i+sample.size)]>0)
    indi=union(d,indi)
  }
}
if (sample.size==2) indi= indi
delta[indi]=1

return(list(models=models_num,delta=delta))
}


###example matrix to calculate PPAs####
freq=cbind(D2=c(0,0,57,0),D3=c(0,0,44,0),D4=c(3,0,45,0))
rownames(freq)=c("a","c","g","t")


#######function to calculate the logP(D|M) given a frequency vector and a mutation rate#####
logPDM=function(z,pstar=.1,mod_num){
  ####z is a frequency vector of length 4, pstar is a mutation rate###
  #p=0
  S=sum(z)
  #comn1=log(pbeta(pstar,S-z[4]+1,z[4]+1,log.p=F))+lbeta(S-z[4]+1,z[4]+1)-log(pstar)
  #comn2=log(1-pbeta(pstar,S-z[4]+1,z[4]+1))+lbeta(S-z[4]+1,z[4]+1)-log(1-pstar)
  comn1=pbeta(pstar,S-z[4]+1,z[4]+1,log.p=T)+lbeta(S-z[4]+1,z[4]+1)-log(pstar)
  comn2=pbeta(pstar,S-z[4]+1,z[4]+1,log.p=T,lower.tail = F)+lbeta(S-z[4]+1,z[4]+1)-log(1-pstar)
  if (mod_num==0) p=(S-z[4])*log(1/3)+comn1
  if (mod_num==1) p=(z[2]+z[3])*log(1/2)+lfactorial(z[1])+lfactorial(z[2]+z[3])-lfactorial(S-z[4]+1)+comn1
  if (mod_num==2) p=(z[1]+z[3])*log(1/2)+lfactorial(z[2])+lfactorial(z[1]+z[3])-lfactorial(S-z[4]+1)+comn1
  if (mod_num==3) p=(z[1]+z[2])*log(1/2)+lfactorial(z[3])+lfactorial(z[1]+z[2])-lfactorial(S-z[4]+1)+comn1
  if (mod_num==4) p=log(2)+lfactorial(z[1])+lfactorial(z[2])+lfactorial(z[3])-lfactorial(S-z[4]+2)+comn1
  if (mod_num==5) p=(S-z[4])*log(1/3)+comn2
  if (mod_num==6) p=(z[2]+z[3])*log(1/2)+lfactorial(z[1])+lfactorial(z[2]+z[3])-lfactorial(S-z[4]+1)+comn2
  if (mod_num==7) p=(z[1]+z[3])*log(1/2)+lfactorial(z[2])+lfactorial(z[1]+z[3])-lfactorial(S-z[4]+1)+comn2
  if (mod_num==8) p=(z[1]+z[2])*log(1/2)+lfactorial(z[3])+lfactorial(z[1]+z[2])-lfactorial(S-z[4]+1)+comn2
  if (mod_num==9) p=log(2)+lfactorial(z[1])+lfactorial(z[2])+lfactorial(z[3])-lfactorial(S-z[4]+2)+comn2
  return(p)
}

#######################sum of logP(D|M) given a model structure vector, a ordered frequency matrix####
sum_logPDM=function(modvec,order.freqmat,pstar=.1){
####modelvec; such as c(7,7,7,0,1,2), c(7,7,7,0,0,0) or c(7,8,7,0,1,0)#####
####the order.freqmat is a ordered freq matrix, i.e., the last row is the mojarity or referenced#####
nsamp=length(modvec)/2
if (nsamp!=dim(order.freqmat)[2])  stop(paste("moder structure and frequency matrix does not match", "\n", ""))
######define 2 data frames###
df.mod =data.frame(mod=modvec[1:nsamp],Index=modvec[-1:-nsamp])
df.freq=data.frame(Index=modvec[-1:-nsamp],freq=t(order.freqmat))
######merge 2 data frames by Index####
df.merged=merge(unique(df.mod),ddply(df.freq,"Index",numcolwise(sum)))
output=0
for (i in 1:dim(df.merged)[1]){
 output[i]=logPDM(z=as.numeric(df.merged[i,3:6]),pstar=pstar,mod_num = df.merged[i,2])
}
return(sum(output))
}



####function to calculate PPAs given a frequency matrix and a reference site for one site####
ppa.gen=function(freqmatrix,pstar=.1,priorPA=.001){
  #freq=freqmatrix
  #if (!any(c("A","C","G","T")==toupper(ref))) ref=NULL
#
 # if (is.null(ref)==T ) {
  #  #####let the majority to be the last row####
   # major.row=which(rank(rowSums(freqmatrix),ties.method = "first")==4)
  #} else {
   # major.row=which(tolower(rownames(freqmatrix))==tolower(ref))
  #}
   # freq[major.row,]=freqmatrix[4,]
    #freq[4,]=freqmatrix[major.row,]
  #  rownames(freq)[major.row]=rownames(freqmatrix)[4]
   # rownames(freq)[4]=rownames(freqmatrix)[major.row]

#####################model structure#######
    x=modmat(dim(freqmatrix)[2])
    model.struc=x$models
    delta=x$delta
    ########vector for p(m_k)####
    pmk=delta*priorPA/sum(delta)+(1-delta)*(1-priorPA)/sum(1-delta)

    ## pmk=c(rep(1-priorPA,10)/10,rep(priorPA,num_interest)/num_interest)

    result.PDM=apply(model.struc,1,FUN = sum_logPDM,order.freqmat=freqmatrix,pstar=pstar)
    PPA=exp(result.PDM)*pmk
    PPA=ifelse(PPA==0,PPA,PPA/sum(PPA))

###############output#####################
    return(list(Model.Top5=model.struc[order(-PPA),][1:5,],PPA.Top5=PPA[order(-PPA)][1:5],PPA_SI=sum(delta*PPA)))

    #return(list(Model=model.struc,
                #log_PDM=result.PDM,
                #PPA=PPA, pmk=pmk,delta=delta))

}





