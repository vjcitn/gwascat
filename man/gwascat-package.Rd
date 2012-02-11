\name{gwascat-package}
\alias{gwascat-package}
\alias{gwascat}
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
Version: \tab 0.0.3\cr
Suggests: \tab \cr
Depends: \tab R (>= 2.14.0), methods, IRanges, GenomicRanges\cr
Imports: \tab \cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
Built: \tab R 2.15.0; ; 2012-02-10 21:08:32 UTC; unix\cr
}

Index:
\preformatted{
gwaswloc-class          Class '"gwaswloc"'
}
Upon attachment, a \code{\link[GenomicRanges]{GRanges-class}} structure
call \code{gwrngs} is formed which can be interrogated by position or
through use of element metadata to learn about catalogued GWAS associations.
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
