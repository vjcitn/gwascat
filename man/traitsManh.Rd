\name{traitsManh}
\alias{traitsManh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
use ggbio facilities to display GWAS results for selected
traits in genomic coordinates
}
\description{
use ggbio facilities to display GWAS results for selected
traits in genomic coordinates
}
\usage{
traitsManh(gwr, selr = GRanges(seqnames = "chr17", IRanges(3e+07, 5e+07)), traits = c("Asthma", "Parkinson's disease", "Height", "Crohn's disease"), truncmlp = 25, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gwr}{
GRanges instance as managed by the gwaswloc container design,
with Disease.Trait and Pvalue\_mlog among elementMetadata columns
}
  \item{selr}{
A GRanges instance to restrict the \code{gwr} for visualization.
Not tested for noncontiguous regions.
}
  \item{traits}{
Character vector of traits to be exhibited; GWAS results
with traits not among these will be labeled ``other''.
}
  \item{truncmlp}{
Maximum value of -log10 p to be displayed; in the raw data
this ranges to the hundreds and can cause bad compression.
}
  \item{\dots}{
not currently used
}
}
\details{
uses a ggbio autoplot
}
\value{
autoplot value
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# do a p-value truncation if you want to reduce compression
traitsManh(gwrngs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ graphics }
