\name{makeCurrentGwascat}
\alias{makeCurrentGwascat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
read NHGRI GWAS catalog table and construct associated GRanges instance
}
\description{
read NHGRI table and construct associated GRanges instance
}
\usage{
makeCurrentGwascat(table.url =
                 "http://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
                 fixNonASCII = FALSE, useHg38seqinfo = TRUE,
                 altSeqinfo, withOnt = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table.url}{
string identifying the .txt file curated at EBI/EMBL
}
  \item{fixNonASCII}{
logical, if TRUE, non-ASCII characters as identified by iconv will be replaced by asterisk
}
  \item{useHg38seqinfo}{logical indicating that metadata on reference
genome comes from GRCh38/hg38}
  \item{altSeqinfo}{Seqinfo instance to be used directly to
document the value produced}
 \item{withOnt}{logical indicating whether 'alternative' (ontology-present,
includes repetition of loci with one:many ontological mapping) or 'full'
(ontology-absent, one record per locus report) version of distributed table}
}
\details{
records for which clear genomic position cannot be determined are dropped from the
ranges instance

an effort is made to use reasonable data types for GRanges metadata, so some qualifying
characters such as (EA) in Risk allele frequency field will simply be omitted during
coercion of contents of that field to numeric.
}
\value{
a GRanges instance
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
# if you have good internet access
newcatr = makeCurrentGwascat()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
