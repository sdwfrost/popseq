% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sfsR.R
\name{sfsR}
\alias{sfsR}
\title{A function to calculate Tajima's D and Fay and WU's H statsitc given a haplotype matrix.}
\usage{
sfsR(hapmatrix)
}
\arguments{
\item{hapmatrix}{A haplotype matrix with ones (ancestral) and twos (derived).}
}
\value{
Tajima's D and Fay and WU's H statsitc.
}
\description{
A function to calculate Tajima's D and Fay and WU's H statsitc given a haplotype matrix.
}
\examples{
exmtx=matrix(sample(2,80,replace=T),ncol=10,nrow=8)
sfsR(exmtx)
sfsR(3-exmtx)


}
\author{
Fei Xiang (\email{xf3087@gmail.com})
}
