#########a function to convert a list of fasta files in to a list by frequencies#######
fastaTolist=function(fas_list,ref){
  #require(ape)
  require(seqinr)
  #dna.list=lapply(fas_list,read.dna,format="fasta",as.character = T,as.matrix = T)
  dna.list=lapply(fas_list,read.alignment,format="fasta")
  dna.list=lapply(dna.list,as.matrix.alignment)


  frequency.list<-lapply(dna.list, function(x) {
    apply(x,2,function(y) table(factor(y, levels = c("a","c","g","t"))))
  })
  ####calculate the consesus based on ref
  cons <- calc_cons(dna.list[[1]], ref)
  #### pstar calculation
  pos=s2n(cons)+1
  ####those consesus position
  pos.cons=cbind(pos,1:dim(frequency.list[[1]])[2])
  freq_sum=Reduce("+", frequency.list)
  pstar=1-sum(freq_sum[pos.cons])/sum(freq_sum)
  ####ordered matrix
  ordered.freq.list=frequency.list
  ordered.freq.list=lapply(ordered.freq.list,function(x) {
    d=x
    d[4,]=x[pos.cons]
    d[pos.cons]=x[4,]
    return(d)
  })
  ##do.call(rbind,ordered.freq.list)
 #### output
  return(list(pstar=pstar,
              freq.mat.ordered=do.call(rbind,ordered.freq.list),
              freq.mat=do.call(rbind,frequency.list),
              cons=cons
  ))

}

