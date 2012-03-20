\name{gwdf_2012_02_02}
\alias{gwdf_2012_02_02}
\alias{gwdf_2012_03_07}
\alias{gwdf_2012_03_20}
\docType{data}
\title{
internal data frame for NHGRI GWAS catalog
}
\description{
convenience container for imported table from NHGRI GWAS catalog
}
\usage{data(gwdf_2012_02_02) # or more recent elements available}
\format{
  A data frame with 7309 observations on the following 34 variables.
  \describe{
    \item{\code{Date Added to Catalog}}{a character vector}
    \item{\code{PUBMEDID}}{a character vector}
    \item{\code{First Author}}{a character vector}
    \item{\code{Date}}{a character vector}
    \item{\code{Journal}}{a character vector}
    \item{\code{Link}}{a character vector}
    \item{\code{Study}}{a character vector}
    \item{\code{Disease/Trait}}{a character vector}
    \item{\code{Initial Sample Size}}{a character vector}
    \item{\code{Replication Sample Size}}{a character vector}
    \item{\code{Region}}{a character vector}
    \item{\code{Chr_id}}{a character vector}
    \item{\code{Chr_pos}}{a character vector}
    \item{\code{Reported Gene(s)}}{a character vector}
    \item{\code{Mapped_gene}}{a character vector}
    \item{\code{Upstream_gene_id}}{a character vector}
    \item{\code{Downstream_gene_id}}{a character vector}
    \item{\code{Snp_gene_ids}}{a character vector}
    \item{\code{Upstream_gene_distance}}{a character vector}
    \item{\code{Downstream_gene_distance}}{a character vector}
    \item{\code{Strongest SNP-Risk Allele}}{a character vector}
    \item{\code{SNPs}}{a character vector}
    \item{\code{Merged}}{a character vector}
    \item{\code{Snp_id_current}}{a character vector}
    \item{\code{Context}}{a character vector}
    \item{\code{Intergenic}}{a character vector}
    \item{\code{Risk Allele Frequency}}{a character vector}
    \item{\code{p-Value}}{a character vector}
    \item{\code{Pvalue_mlog}}{a character vector}
    \item{\code{p-Value (text)}}{a character vector}
    \item{\code{OR or beta}}{a character vector}
    \item{\code{95\% CI (text)}}{a character vector}
    \item{\code{Platform..}}{a character vector}
    \item{\code{CNV}}{a character vector}
  }
}
\source{
\url{http://www.genome.gov/gwastudies}
}
\note{
The .onAttach function specifies which data frame is transformed
to GRanges.
}
%\references{
%}
\examples{
\dontrun{
data(gwdf_2012_02_02)
data(gwdf_2012_03_07)
}
}
\keyword{datasets}
