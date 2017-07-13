#' An efficient function to identify sites of interestes from a list of longitudinal sequence file.
#'
#' @param file_list A list of bam files.
#' @param ref A referece sequence.
#' @param priorPA Prior assocation to possible models
#' @param region An optional region from the BAM file (as a \code{GRanges} object)
#' @return PPA association for each site.
#' @examples
#'
#'
#'
#'p4305.d3=system.file("extdata/4305-d3.bam",package = "popseq")
#'p4305.d4=system.file("extdata/4305-d4.bam",package = "popseq")
#'p4305.d5=system.file("extdata/4305-d5.bam",package = "popseq")
#'ref.seg1=system.file("extdata/seg1.fas",package = "popseq")
#'ref.seg2=system.file("extdata/seg2.fas",package = "popseq")
#'
#'
#'pigbam=list(p4305.d3,p4305.d4,p4305.d5)
#'
#'
#'res1=screen_bam_eff(file_list=pigbam,ref=ref.seg1,region= GRanges("GQ166656_1",IRanges(1,2280)),priorPA=.001)
#'
#'res2=screen_bam_eff(file_list=pigbam,ref=ref.seg2,priorPA=.001,region = GRanges("GQ166655_2",IRanges(1,2274)))
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
screen_bam_eff=function(file_list,ref=NULL,region,priorPA=.001){
  #########seq from a list of fasta files######
  seq.fas=bamTolist(bam_list = file_list,ref=ref,region=region)
  pstar=seq.fas$pstar
  freq.mat=seq.fas$freq.mat.ordered
  cons=seq.fas$cons
  nsamp=dim(freq.mat)[1]/4
  z=modmat(nsamp)

  ######calculate the ppas########
  site.PPA_SI=0
  site.PPA_SI=apply(freq.mat,2, function(x) {
    M=matrix(x,nrow=4)
    result=result=ppa.gen_eff(freqmatrix =  M,pstar=pstar,priorPA=priorPA,models =z$models,delta=z$delta )
    return(result$PPA_SI)
  })
  names(site.PPA_SI)=NULL
  site.and.PPA_SI=cbind(site=1:dim(freq.mat)[2],PPA_SI=site.PPA_SI)
  output=list(pstar=pstar,basedist=seq.fas$freq.mat,
              cons=cons,site.PPA_SI=t(site.and.PPA_SI),
              site.interest=site.and.PPA_SI[site.and.PPA_SI[,2]>.5,])
  return(output)

}
