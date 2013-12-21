\name{gwascat-package}
\alias{gwascat-package}
\alias{gwascat}
\alias{g17SM} 
\alias{gg17N} 
\alias{gw6.rs_17} 
\alias{low17}
\alias{rules_6.0_1kg_17}
\alias{DataFrame}
\alias{gwrngs}
\alias{impute.snps}
\alias{SnpMatrix-class}
\docType{package}
\title{
representing and modeling data in the NHGRI GWAS catalog
}
\description{
representing and modeling data in the NHGRI GWAS catalog,
using GRanges and allied infrastructure
}
\details{
\tabular{ll}{
Package: \tab gwascat\cr
Version: \tab 1.7.3\cr
Suggests: \tab \cr
Depends: \tab R (>= 3.0.0), methods, IRanges, GenomicRanges\cr
Imports: \tab \cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
gwaswloc-class          Class '"gwaswloc"'
}
Upon attachment, a \code{\link[GenomicRanges]{GRanges-class}} structure
call \code{gwrngs} is formed which can be interrogated by position or
through use of element metadata to learn about catalogued GWAS associations.

The data objects

'g17SM' 'gg17N' 'gw6.rs_17' 'low17' 'rules_6.0_1kg_17'
'gwrngs'

are described in vignettes.

The DataFrame function is imported from IRanges.

The \code{\link[snpStats]{SnpMatrix-class}} is used to represent data
related to rule-based imputation, using the \code{\link[snpStats]{impute.snps}}
function.
}
\author{
VJ Carey <stvjc@channing.harvard.edu>

Maintainer: VJ Carey <stvjc@channing.harvard.edu>
}
\references{
\url{http://www.genome.gov/gwastudies/}.

Partial support from the Computational Biology Group at
Genentech, Inc.

}
\keyword{ package }
%\seealso{
%}
\examples{
\dontrun{
 gwrngs
 }
}
