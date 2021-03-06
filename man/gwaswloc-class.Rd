\name{gwaswloc-class}
\Rdversion{1.1}
\docType{class}
\alias{gwaswloc-class}
\alias{[,gwaswloc,ANY,ANY,ANY-method}
\alias{[,gwaswloc,ANY-method}
\alias{[,gwaswloc-method}
\alias{show,gwaswloc-method}
\alias{getRsids,gwaswloc-method}
\alias{getTraits,gwaswloc-method}
\alias{subsetByChromosome,gwaswloc-method}
\alias{subsetByTraits,gwaswloc-method}
\alias{getTraits}
\alias{getRsids}
\alias{subsetByChromosome}
\alias{subsetByTraits}


\title{Class \code{"gwaswloc"}}
\description{
A container for GRanges instances representing information in the NHGRI GWAS catalog.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gwaswloc", ...)}.
Any GRanges instance can be supplied.
}
\section{Slots}{
  \describe{
    \item{\code{extractDate}:}{character set manually in .onAttach code to
indicate date of retrieval of base table}
    \item{\code{seqnames}:}{Object of class \code{"Rle"} typically
representing chromosome numbers of loci associated with specific traits }
    \item{\code{ranges}:}{Object of class \code{"IRanges"} genomic
coordinates for locus }
    \item{\code{strand}:}{Object of class \code{"Rle"} identifier of
chromosome strand }
    \item{\code{elementMetadata}:}{Object of class \code{"DataFrame"} general
\code{\link{DataFrame-class}} instance providing attributes
for the locus-trait association }
    \item{\code{seqinfo}:}{Object of class \code{"Seqinfo"}  }
    \item{\code{metadata}:}{Object of class \code{"list"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{GRanges}"}, directly.
Class \code{"\linkS4class{GenomicRanges}"}, by class "GRanges", distance 2.
Class \code{"\linkS4class{Vector}"}, by class "GRanges", distance 3.
Class \code{"\linkS4class{GenomicRanges_OR_missing}"}, by class "GRanges", distance 3.
Class \code{"\linkS4class{GenomicRanges_OR_GRangesList}"}, by class "GRanges", distance 3.
Class \code{"\linkS4class{Annotated}"}, by class "GRanges", distance 4.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "gwaswloc")}: a character argument to the bracket
will be assumed to be a dbSNP identifier for a SNP locus, and records corresponding
to this SNP are extracted; numeric indexes are supported as for
\code{\link{GRanges-class}} instances.  }
    \item{getRsids}{\code{signature(x = "gwaswloc")}: extract all dbSNP identifiers
as a character vector}
    \item{getTraits}{\code{signature(x = "gwaswloc")}: extract all traits (NHGRI
term 'Disease/Trait') as a
character vector}
    \item{subsetByChromosome}{\code{signature(x = "gwaswloc")}: select records by chromosome, a vector of chromosomes may be supplied }
    \item{subsetByTraits}{\code{signature(x = "gwaswloc")}: select all records corresponding
to a given vector of traits }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
\url{http://www.ebi.ac.uk/gwas/}
}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
In gwascat package 1.9.6 and earlier, 
the globally accessible \code{gwaswloc} instance
\code{gwrngs} was created upon attachment.  This is no 
longer the case.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("gwaswloc")
}
\keyword{classes}
