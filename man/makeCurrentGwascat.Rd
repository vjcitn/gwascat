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
                 fixNonASCII = TRUE, genome="GRCh38",
                 withOnt = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table.url}{
string identifying the .txt file curated at EBI/EMBL
}
  \item{fixNonASCII}{
logical, if TRUE, non-ASCII characters as identified by iconv will be replaced by asterisk
}
 \item{genome}{character string: 'GRCh38' is default and yields
current image as provided by EMBL/EBI; 'GRCh37' yields a realtime liftOver
to hg19 coordinates, via AnnotationHub storage of the chain files.
Any other value yields an error.}
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
# if you have good internet access
  if (interactive()) {
     newcatr = makeCurrentGwascat()
     newcatr
     }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
