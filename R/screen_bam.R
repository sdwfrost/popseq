#' A function to identify sites of interestes from a list of longitudinal sequence file.
#'
#' @param file_list A list of bam files.
#' @param ref A referece sequence.
#' @param priorPA Prior assocation to possible models
#' @param region An optional region from the BAM file (as a \code{GRanges} object)
#' @return PPA association for each site.
#' @examples
#' pig2.1=system.file("extdata/pig2.1",package = "popseq")
#' pig2.2=system.file("extdata/pig2.2",package = "popseq")
#' pig2.3=system.file("extdata/pig2.3",package = "popseq")
#'
#' ref=system.file("extdata/regions.fas",package = "popseq")
#'
#' pig2=list(pig2.1,pig2.2,pig2.3)
#'
#' screen_select(pig2,ref=ref)
#'
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#'
#' @export


#########an overall function to deal with a list of bam files on selection#######
screen_bam=function(file_list,ref=NULL,region,priorPA=.001){
  #########seq from a list of fasta files######
  seq.fas=bamTolist(bam_list = file_list,ref=ref,region=region)
  pstar=seq.fas$pstar
  freq.mat=seq.fas$freq.mat.ordered
  cons=seq.fas$cons
  ######calculate the ppas########
  site.PPA_SI=0
  site.PPA_SI=apply(freq.mat,2, function(x) {
    M=matrix(x,nrow=4)
    result=ppa.gen(freqmatrix =  M,pstar=pstar,priorPA=priorPA)
    return(result$PPA_SI)
  })
  names(site.PPA_SI)=NULL
  site.and.PPA_SI=cbind(site=1:dim(freq.mat)[2],PPA_SI=site.PPA_SI)
  output=list(pstar=pstar,basedist=seq.fas$freq.mat,
              cons=cons,site.PPA_SI=t(site.and.PPA_SI),
              site.interest=site.and.PPA_SI[site.and.PPA_SI[,2]>.5,])
  return(output)

}
