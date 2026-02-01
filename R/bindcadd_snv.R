#' bind CADD scores of Kircher et al. 2014 to a GRanges instance
#' 
#' bind CADD scores of Kircher et al. 2014 to a GRanges instance; by default
#' will use HTTP access at UW
#' 
#' joins CADD fields at addresses that match query; the CADD fields for query
#' ranges that are not matched are set to NA
#' 
#' @param gr query ranges to which CADD scores should be bound
#' @param fn path to Tabix-indexed bgzipped TSV of CADD as distributed at
#' krishna.gs.washington.edu on 1 April 2014
#' @return GRanges instance with additional fields as obtained in the CADD
#' resource
#' @note This software developed in part with support from Genentech, Inc.
#' @author VJ Carey <stvjc@@channing.harvard.edu>
#' @references M Kircher, DM Witten, P Jain, BJ O'Roak, GM Cooper, J Shendure,
#' A general framework for estimating the relative pathogenicity of human
#' genetic variants, Nature Genetics Feb 2014, PMID 24487276
#' @keywords models
#' @examples
#' 
#'  \dontrun{
#'   data(ebicat_2020_04_30)
#'   g2 = as(ebicat_2020_04_30, "GRanges")
#'  # would need to lift over here 
#'   bindcadd_snv( g2[which(seqnames(g2)=="chr2")][1:20] )
#'   }
#' 
#' @export bindcadd_snv
bindcadd_snv = function(gr, fn="http://krishna.gs.washington.edu/download/CADD/v1.0/1000G.tsv.gz") {
#
# import and export SNVs from CADD
#
#  subroutine to set up
#  mcols for result
#
#  want DataFrame of NA to be populated only with attributes of matched
#  ranges
padDF = function(n, targDF) {
    nna = rep(NA, n)
    nc = ncol(targDF)
    cl = sapply(targDF, class)
    if (any(cl=="factor")) cl[cl=="factor"] = "character"
    exemplars = lapply(lapply(cl, function(x) as(NA, x)), rep, n)
    ini = DataFrame(exemplars)
    names(ini) = names(targDF)
    ini
    }
# end subroutine
#
 tf = TabixFile(fn)
 open(tf)
 on.exit(close(tf))
 seqlevels(gr) = gsub("chr", "", seqlevels(gr))
 cur = import(tf, which=gr)
 ncc = function(x) nchar(as.character(x))
 issn = which(ncc(cur$V3)==1 & ncc(cur$V4)==1)  # restrict to SNV
 names(values(cur)) = c("Ref", "Alt", "CScore", "PHRED")
 newdf = padDF(length(gr), mcols(cur))
 fo = findOverlaps(gr, cur)
  for (i in 1:ncol(newdf)) {
        f = force
        if (class(mcols(cur)[,i]) == "factor") f = function(x) as.character(x)
        newdf[queryHits(fo), i] = f(mcols(cur)[subjectHits(fo), i])
        }
    mcols(gr) = cbind(mcols(gr), newdf)
 gr
}

