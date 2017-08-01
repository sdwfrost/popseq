#########a function to convert a list of bam files in to a list by frequencies#######
bamTolist=function(bam_list,ref,region){
  require(seqinr)
  #######read the full ref and select the region matching segment#######
  whole.ref=read.alignment(ref,format="fasta")
  piece=match(region@seqnames@values,whole.ref$nam)
  dna.ref=strsplit(whole.ref$seq[[piece]],"")
  dna.ref=matrix(unlist(dna.ref),nrow=1)
  #dna.ref=tolower(as.matrix.alignment(read.alignment(ref,format="fasta")))
  segment.bam <- ScanBamParam(which=region)


  bamdna.list=lapply(bam_list, function(x) {
    xx=BamFile(x)
    piup=pileup(xx,scanBamParam = segment.bam)
    freq=pileupFreq(piup)
    t(freq[,3:6])
  })

 ####produce the consesus#####

  #####consesus from the first bam file###

  cons.bam1=apply(bamdna.list[[1]],2,function(x) rownames(bamdna.list[[1]])[order(-x)][1])
  cons.bam1=tolower(cons.bam1)
###########some warnings needed here##########
  cons=apply(dna.ref,2,function(x) {
    if (x %in% c("a","c","g","t")){
      return(x)
    } else return (NA)
  })

  names(cons)=NULL
  cons[is.na(cons)]=cons.bam1[is.na(cons)]
  ####same process then on as shown in fastaTolist#####
  #### pstar calculation
  pos=s2n(cons)+1
  ####those consesus position
  pos.cons=cbind(pos,1:dim(bamdna.list[[1]])[2])
  freq_sum=Reduce("+", bamdna.list)
  pstar=1-sum(freq_sum[pos.cons])/sum(freq_sum)

  ####ordered matrix
  ordered.freq.list=bamdna.list
  ordered.freq.list=lapply(ordered.freq.list,function(x) {
    d=x
    d[4,]=x[pos.cons]
    d[pos.cons]=x[4,]
    return(d)
  })

  #### output
  return(list(pstar=pstar,
              freq.mat.ordered=do.call(rbind,ordered.freq.list),
              freq.mat=do.call(rbind,bamdna.list),
              cons=cons
  ))


}



