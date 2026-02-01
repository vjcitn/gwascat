#' container for gwas hit data and GRanges for addresses
#' @import methods
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
setClass("gwaswloc", representation(extractDate="character"), 
   contains="GRanges")

setMethod("show", "gwaswloc", function(object) {
 cat("gwasloc instance with", length(object), "records and", ncol(mcols(object)),
  "attributes per record.\n")
 cat("Extracted: ", object@extractDate, "\n")
 if ("badpos" %in% names(S4Vectors::metadata(object))) {
   cat("metadata()$badpos includes records for which no unique locus was given.\n")
   }
 if (length(object)==0) {
     cat("no records selected.\n")
     return(invisible(NULL))
     }
 cat("Genome: ", unique(genome(object)), "\n")
 cat("Excerpt:\n")
 nrec = min(5, length(object))
 availcols = colnames(mcols(object))
 majorfields = c("DISEASE/TRAIT", "SNPS", "P-VALUE")
 if (all(majorfields %in% availcols) & length(availcols) > 5)
   show(as(object, "GRanges")[1:nrec, majorfields])
 else show(as(object, "GRanges")[1:nrec,])
})

#
# intention here is to allow row subscripting by rs number without any
# casting ... and to allow all other methods to proceed as in GRanges API
#
#' extractor for gwaswloc
#' @param x gwaswloc
#' @param i index
#' @param j index
#' @param \dots addtl indices
#' @param drop logical(1)
#' @export
setMethod("[", "gwaswloc", function(x, i, j, ..., drop=FALSE) {
 if (missing(drop)) drop <- FALSE
# if (!missing(j)) stop("no column subscripting on gwaswloc")
 if (missing(i)) i = 1:length(x)
 if (!is(i, "numeric")) {
  rsids = getRsids(x)
  i = match(i, rsids, nomatch=0)
 }
 if (isTRUE(length(i)==0)) stop("index has length 0")
 if (isTRUE(i==0)) return(new("gwaswloc"))
 callNextMethod()
})
 

#' generic snp name retrieval
#' @param x gwaswloc
#' @export
setGeneric("getRsids", function(x)standardGeneric("getRsids"))
#' specific snp name retrieval
#' @param x gwaswloc
#' @export
setMethod("getRsids", "gwaswloc", function(x)
 mcols(x)$SNPS)

#' generic trait retrieval
#' @param x gwaswloc
#' @export
setGeneric("getTraits", function(x)standardGeneric("getTraits"))
#' specific trait retrieval
#' @param x gwaswloc
#' @export
setMethod("getTraits", "gwaswloc", function(x)
 mcols(x)[["DISEASE/TRAIT"]])


#' generic trait subsetting
#' @param x gwaswloc
#' @param ch character vector of chromosomes
#' @export
setGeneric("subsetByChromosome", function(x, ch)standardGeneric("subsetByChromosome"))
#' specific trait subsetting
#' @param x gwaswloc
#' @param ch character vector of chromosomes
#' @export
setMethod("subsetByChromosome", "gwaswloc", function(x, ch) {
 x[ which(as.character(seqnames(x)) %in% ch) ]
})

#' generic trait subsetting
#' @param x gwaswloc
#' @param tr character vector of traits
#' @export
setGeneric("subsetByTraits", function(x, tr)standardGeneric("subsetByTraits"))
#' specific trait subsetting
#' @param x gwaswloc
#' @param tr character vector of traits
#' @export
setMethod("subsetByTraits", "gwaswloc", function(x, tr) {
 x[ which(getTraits(x) %in% tr) ]
})

