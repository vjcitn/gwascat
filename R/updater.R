
makeCurrentGwascat = function(table.url="http://www.genome.gov/admin/gwascatalog.txt", fixNonASCII=TRUE, useHg38seqinfo = TRUE, altSeqinfo) {
 tab = read.delim(url(table.url), sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
 if (missing(altSeqinfo) & !useHg38seqinfo) stop("need an altSeqinfo when  useHg38seqinfo is FALSE")
 if (fixNonASCII) tab = fixNonASCII(tab)
 cur = gwdf2GRanges(tab, extractDate=as.character(Sys.Date()))
 seqlevelsStyle(cur) = "NCBI"
 cursn = seqlevels(cur)
 if (useHg38seqinfo) {
    data(si.hs.38)  
    seqinfo(cur) = si.hs.38[cursn]
    }
 else seqinfo(cur) = altSeqinfo[cursn]
 cur
}
