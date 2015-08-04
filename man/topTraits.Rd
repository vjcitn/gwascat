\name{topTraits}

\alias{topTraits}
\alias{chklocs}
\alias{locs4trait}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
operations on GWAS catalog
}
\description{
operations on GWAS catalog
}
\usage{
topTraits (gwwl, n=10, tag="DISEASE.TRAIT") 

locs4trait(gwwl, trait, tag="DISEASE.TRAIT") 

chklocs(chrtag="20", gwwl)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gwwl}{
instance of \code{\linkS4class{gwaswloc}}
}
 \item{n}{numeric, number of traits to report}
 \item{tag}{character, name of field to be used for trait enumeration}
 \item{trait}{character, trait to use for filtering}
 \item{chrtag}{character, chromosome identifier}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{topTraits} returns a character vector of most frequently
occurring traits in the database

\code{locs4trait} returns a \code{\linkS4class{gwaswloc}}
object with records defining associations to the specified trait

\code{chklocs} returns a logical that is TRUE when the
asserted locations of SNP in the GWAS catalog agree
with the locations given in the dbSNP package 
SNPlocs.Hsapiens.dbSNP.20110815
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
data(ebicat38)
topTraits(ebicat38)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
