\name{bindcadd_snv}
\alias{bindcadd_snv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bind CADD scores of Kircher et al. 2014 to a GRanges instance
}
\description{
bind CADD scores of Kircher et al. 2014 to a GRanges instance; by
default will use HTTP access at UW
}
\usage{
bindcadd_snv(gr, fn = "http://krishna.gs.washington.edu/download/CADD/v1.0/1000G.tsv.gz")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gr}{
query ranges to which CADD scores should be bound
}
  \item{fn}{
path to Tabix-indexed bgzipped TSV of CADD as distributed at krishna.gs.washington.edu on
1 April 2014
}
}
\details{
joins CADD fields at addresses that match query; the CADD fields for
query ranges that are not
matched are set to NA
}
\value{
GRanges instance with additional fields as obtained in the CADD resource
}
\references{
M Kircher, DM Witten, P Jain, BJ O'Roak, GM Cooper, J Shendure,
A general framework for estimating the relative pathogenicity of human genetic variants,
Nature Genetics Feb 2014, PMID 24487276
}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
This software developed in part
with support from Genentech, Inc.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
 \dontrun{
 # requires internet access
 data(gwrngs19)
 g2 = as(gwrngs19, "GRanges")
 bindcadd_snv( g2[which(seqnames(g2)=="chr2")][1:20] )
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
