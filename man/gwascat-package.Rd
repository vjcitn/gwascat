\name{gwascat-package}
\alias{gwascat-package}
\alias{gwascat}
\alias{g17SM} 
\alias{gg17N} 
\alias{gw6.rs_17} 
\alias{low17}
\alias{rules_6.0_1kg_17}
\alias{DataFrame}
\alias{gwrngs19}
\alias{gwrngs38}
\alias{ebicat38}
\alias{ebicat37}
\alias{impute.snps}
\alias{SnpMatrix-class}
\alias{si.hs.38}
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

The GWAS catalog management has migrated to EMBL/EBI.  Use 
data(ebicat38) for an image dated 3 August 2015.  Use makeCurrentGwascat()
to get a more recent image.  Use data(ebicat37) for a GRCh37 (or hg19)
liftOver result.  Use data(ebicat37UCSC) for an image with
\code{hg19} as genome tag and UCSC seqnames.

The data objects

'g17SM' 'gg17N' 'gw6.rs_17' 'low17' 'rules_6.0_1kg_17'
'gwrngs'

are described in vignettes.

The DataFrame function is imported from IRanges.

The \code{\link[snpStats]{SnpMatrix-class}} is used to represent data
related to rule-based imputation, using the \code{\link[snpStats]{impute.snps}}
function.

si.hs.38 is a \code{\link[GenomeInfoDb]{Seqinfo-class}} instance for hg38.
}
\author{
VJ Carey <stvjc@channing.harvard.edu>

Maintainer: VJ Carey <stvjc@channing.harvard.edu>
}
\references{
\url{http://www.ebi.ac.uk/gwas/}

Partial support from the Computational Biology Group at
Genentech, Inc.

}
\keyword{ package }
\examples{
 data(ebicat38)
 ebicat38
}
