\name{ldtagr}
\alias{ldtagr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
expand a list of variants by including those in a VCF with LD exceeding
some threshold
}
\description{
expand a list of variants by including those in a VCF with LD exceeding
some threshold
}
\usage{
ldtagr(snprng, tf, samples, genome = "hg19", lbmaf = 0.05, lbR2 = 0.8, radius = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snprng}{
a named GRanges for a single SNP.  The name must correspond to
the name that will be assigned by \code{\link[VariantAnnotation]{genotypeToSnpMatrix}} to the corresponding column of a SnpMatrix.
}
  \item{tf}{
TabixFile instance pointing to a bgzipped tabix-indexed VCF file
}
  \item{samples}{
a vector of sample identifiers, if excluded, all samples used
}
  \item{genome}{
tag like 'hg19'
}
  \item{lbmaf}{
lower bound on variant MAF to allow consideration
}
  \item{lbR2}{
lower bound on R squared for regarding SNP to be incorporated
}
  \item{radius}{
radius of search in bp around the input range
}
}
\details{
uses snpStats ld()
}
\value{
a GRanges with names corresponding to 'new' variants
and mcols fields 'paramRangeID' (base variant input)
and 'R2'
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey
}
\note{
slow but safe approach.  probably a matrix method could
be substituted using the nice sparse approach already in snpStats
}


%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
require(GenomicRanges)
cand = GRanges("1", IRanges(113038694, width=1))
names(cand) = "rs883593"
require(VariantAnnotation)
expath = dir(system.file("vcf", package="gwascat"), patt=".*exon.*gz$", full=TRUE)
tf = TabixFile(expath)
ldtagr( cand, tf, lbR2 = .8)
# should do with 1000 genomes in S3 bucket and gwascat
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
