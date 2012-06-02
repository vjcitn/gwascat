\name{obo2graphNEL}
\alias{obo2graphNEL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
convert a typical OBO text file to a graphNEL instance (using Term elements)
}
\description{
convert a typical OBO text file to a graphNEL instance (using Term elements)
}
\usage{
obo2graphNEL(obo, kill = "\\\\[Typedef\\\\]")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obo}{
string naming a file in OBO format
}
  \item{kill}{
entity types to be excluded from processing -- probably this
should be in a 'keep' form, but for now this works.
}
}
\details{
Very rudimentary list and grep operations are used to
retain sufficient information to map the DAG to a graphNEL,
using formal term identifiers as node names and 'is-a' relationships
as edges, and term names and other metadata are
assigned to nodeData components.
}
\value{
a graphNEL instance
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
For use with human disease ontology,
\url{http://www.obofoundry.org/cgi-bin/detail.cgi?id=disease_ontology}
}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
The OBO for Human Disease ontology is serialized as text with this
package.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
