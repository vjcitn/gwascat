\name{gwcex2gviz}
\alias{gwcex2gviz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare salient components of GWAS catalog for rendering with Gviz
}
\description{
Prepare salient components of GWAS catalog for rendering with Gviz
}
\usage{
gwcex2gviz(contextGR = 
  GRanges(seqnames = "chr17", IRanges(start = 37500000, width = 1e+06)), 
  txrefpk = "TxDb.Hsapiens.UCSC.hg19.knownGene", 
  genesympk = "org.Hs.eg.db",
  genome = "hg19", plot.it = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contextGR}{
A GRanges instance delimiting the visualization in genomic coordinates
}
  \item{txrefpk}{
a TxDb package, typically
}
  \item{genesympk}{
string naming annotationDbi .db package
}
  \item{genome}{
character tag like 'hg19'
}
  \item{plot.it}{
logical, if FALSE, just return list
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%%  ~Describe the value returned
%%%  If it is a LIST, use
%%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
args(gwcex2gviz)
gwcex2gviz()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphics }