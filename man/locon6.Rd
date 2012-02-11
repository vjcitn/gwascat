\name{locon6}
\alias{locon6}
\docType{data}
\title{
location information for 10000 SNPs probed on Affy GW 6.0
}
\description{
location information for 10000 SNPs probed on Affy GW 6.0
}
\usage{data(locon6)}
\format{
  A data frame with 10000 observations on the following 3 variables.
  \describe{
    \item{\code{dbsnp_rs_id}}{a character vector}
    \item{\code{chrom}}{a character vector}
    \item{\code{physical_pos}}{a numeric vector}
  }
}
\details{
extracted from pd.genomewidesnp.6 v 1.4.0; for demonstration
purposes
}
%\source{
%%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(locon6)
str(locon6) 
}
\keyword{datasets}
