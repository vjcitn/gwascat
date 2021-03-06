\name{gwdf_2012_02_02}
\alias{gwdf_2014_09_08}
\docType{data}
\title{
internal data frame for NHGRI GWAS catalog
}
\description{
convenience container for imported table from NHGRI GWAS catalog
}
\usage{data("gwdf_2014_09_08")}
\format{
  A data frame with 17832 observations on the following 34 variables.
  \describe{
    \item{\samp{Date Added to Catalog}}{a character vector}
    \item{\code{PUBMEDID}}{a numeric vector}
    \item{\samp{First Author}}{a character vector}
    \item{\code{Date}}{a character vector}
    \item{\code{Journal}}{a character vector}
    \item{\code{Link}}{a character vector}
    \item{\code{Study}}{a character vector}
    \item{\samp{Disease/Trait}}{a character vector}
    \item{\samp{Initial Sample Size}}{a character vector}
    \item{\samp{Replication Sample Size}}{a character vector}
    \item{\code{Region}}{a character vector}
    \item{\code{Chr_id}}{a character vector}
    \item{\code{Chr_pos}}{a character vector}
    \item{\samp{Reported Gene(s)}}{a character vector}
    \item{\code{Mapped_gene}}{a character vector}
    \item{\code{Upstream_gene_id}}{a character vector}
    \item{\code{Downstream_gene_id}}{a character vector}
    \item{\code{Snp_gene_ids}}{a character vector}
    \item{\code{Upstream_gene_distance}}{a character vector}
    \item{\code{Downstream_gene_distance}}{a character vector}
    \item{\samp{Strongest SNP-Risk Allele}}{a character vector}
    \item{\code{SNPs}}{a character vector}
    \item{\code{Merged}}{a character vector}
    \item{\code{Snp_id_current}}{a character vector}
    \item{\code{Context}}{a character vector}
    \item{\code{Intergenic}}{a character vector}
    \item{\samp{Risk Allele Frequency}}{a character vector}
    \item{\samp{p-Value}}{a character vector}
    \item{\code{Pvalue_mlog}}{a character vector}
    \item{\samp{p-Value (text)}}{a character vector}
    \item{\samp{OR or beta}}{a character vector}
    \item{\samp{95\% CI (text)}}{a character vector}
    \item{\samp{Platform [SNPs passing QC]}}{a character vector}
    \item{\code{CNV}}{a character vector}
  }
}
\source{
\url{http://www.ebi.ac.uk/gwas/}
}
\note{
In versions prior to 1.9.6,
The .onAttach function specifies which data frame is transformed
to GRanges.  This is now managed manually.
}
%\references{
%}
\examples{
\dontrun{
data(gwdf_2014_09_08)
# try gwascat:::gwdf2GRanges on this data.frame
}
}
\keyword{datasets}
