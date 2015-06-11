
makeCurrentGwascat.legacy = function(table.url="http://www.genome.gov/admin/gwascatalog.txt", fixNonASCII=TRUE, useHg38seqinfo = TRUE, altSeqinfo) {
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

makeCurrentGwascat = function(table.url=
  "http://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
   fixNonASCII=FALSE, useHg38seqinfo = TRUE, altSeqinfo, withOnt=FALSE) {
 suppressWarnings({
 if (!withOnt) table.url = sub("alternative", "full", table.url)
 message(paste0("running read.delim on ", table.url, "..."))
 tab = read.delim(url(table.url), sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
 })
 message(paste0("formatting gwaswloc instance..."))
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
 metadata(cur) = list(
    date.created = date(),
    creation = match.call(),
    sessInfo.creation = sessionInfo()
    )
 message("done.")
 cur
}
