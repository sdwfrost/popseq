#' An efficient function to identify sites of interestes from a list of longitudinal sequence files with parallel computing######
#'
#' @param file_list A list of bam files.
#' @param ref A referece sequence.
#' @param priorPA Prior assocation to possible models.
#' @param region An optional region from the BAM file (as a \code{GRanges} object).
#' @param core Number of cores for parallel computation.
#' @return PPA association for each site.
#' @examples
#'
#'
#'
#'p4305.d3=system.file("extdata/4305-d3.bam",package = "popseq")
#'p4305.d4=system.file("extdata/4305-d4.bam",package = "popseq")
#'p4305.d5=system.file("extdata/4305-d5.bam",package = "popseq")
#'ref=system.file("extdata/A_England_195_2009.fa",package = "popseq")
#'
#'
#'
#'pigbam=list(p4305.d3,p4305.d4,p4305.d5)
#'
#'
#'res1=screen_bam_multi(file_list=pigbam,ref=ref,region= GRanges("GQ166656_1",IRanges(1,2280)),priorPA=.001,core=8)
#'
#'res2=screen_bam_multi(file_list=pigbam,ref=ref,region = GRanges("GQ166655_2",IRanges(1,2274)),priorPA=.001,core=8)
#'
#'
#'
#'
#'
#'
#'
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#'
#' @export


#########an overall function to deal with a list of bam files on selection#######
screen_bam_multi=function(file_list,ref=NULL,region,priorPA=.001,core=1){
  #########seq from a list of fasta files######
  seq.fas=bamTolist(bam_list = file_list,ref=ref,region=region)
  pstar=seq.fas$pstar
  freq.mat=seq.fas$freq.mat.ordered
  cons=seq.fas$cons
  nsamp=dim(freq.mat)[1]/4
  z=modmat(nsamp)

  ######calculate the ppas########
  site.PPA_SI=list()
  freq.mat.list=as.list(data.frame(t(freq.mat)))
  site.PPA_SI=mclapply(freq.mat.list, function(x) {
    M=matrix(x,nrow=4)
    result=ppa.gen_eff(freqmatrix =  M,pstar=pstar,priorPA=priorPA,models =z$models,delta=z$delta )
    return(result$PPA_SI)
  },mc.cores = core)

  site.PPA_SI.vector=unlist(site.PPA_SI)
  names(site.PPA_SI.vector)=NULL
  site.and.PPA_SI=cbind(site=1:dim(freq.mat)[2],PPA_SI=site.PPA_SI.vector)
  output=list(pstar=pstar,basedist=seq.fas$freq.mat,
              cons=cons,site.PPA_SI=t(site.and.PPA_SI),
              site.interest=site.and.PPA_SI[site.and.PPA_SI[,2]>.5,])
  return(output)

}
