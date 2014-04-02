
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

