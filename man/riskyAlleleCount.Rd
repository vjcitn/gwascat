\name{riskyAlleleCount}
\alias{riskyAlleleCount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
given a matrix of subjects x SNP calls, count number of risky alleles
}
\description{
given a matrix of subjects x SNP calls, count number of risky alleles for
various conditions, relative
to NHGRI GWAS catalog
}
\usage{
riskyAlleleCount(callmat, matIsAB = TRUE, chr, 
  gwwl = gwrngs, snpap = "SNPlocs.Hsapiens.dbSNP.20110815", 
  gencode = c("A/A", "A/B", "B/B"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{callmat}{
matrix with subjects as rows, SNPs as columns; entries can be generic
A/A, A/B, B/B, or specific nucleotide calls
}
  \item{matIsAB}{
logical, FALSE if nucleotide codes are present, TRUE if 
generic call codes are present; in the latter case, gwascat:::ABmat2nuc
will be run 
}
  \item{chr}{
code for chromosome, should work with the SNP annotation getSNPlocs
function, so likely "ch[nn]"
}
  \item{gwwl}{
an instance of \code{\linkS4class{gwaswloc}}
}
  \item{snpap}{
name of a Bioconductor SNPlocs.Hsapiens.dbSNP.* package
}
  \item{gencode}{
codes used for generic SNP call
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
matrix with rows corresponding to subjects , columns corresponding
to SNP
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(gg17N)
h17 = riskyAlleleCount(gg17N, matIsAB=FALSE, chr="chr17")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
