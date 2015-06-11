
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
   fixNonASCII=FALSE, genome="GRCh38", withOnt=FALSE) {
 stopifnot(genome %in% c("GRCh37", "GRCh38"))
 suppressWarnings({
 if (!withOnt) table.url = sub("alternative", "full", table.url)
 message(paste0("running read.delim on ", table.url, "..."))
 tab = read.delim(url(table.url), sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
 })
 message(paste0("formatting gwaswloc instance..."))
 if (fixNonASCII) tab = fixNonASCII(tab)
 cur = gwdf2GRanges(tab, extractDate=as.character(Sys.Date()))
 seqlevelsStyle(cur) = "NCBI"
 cursn = seqlevels(cur)
 data(si.hs.38)  
 seqinfo(cur) = si.hs.38[cursn]
 if (genome == "GRCh37") cur = lo38to19(cur)
 metadata(cur) = list(
    date.created = date(),
    creation = match.call(),
    sessInfo.creation = sessionInfo()
    )
 message("done.")
 cur
}

lo38to19 = function(gwwl) {
 message("starting liftover from GRCh38 to GRCh37")
 stopifnot(genome(gwwl)[1] == "GRCh38")
 ah = AnnotationHub()
 ii = query(ah, "UCSC liftOver chain file from hg38 to hg19")
 ch = ah[[names(ii)]]
 seqlevelsStyle(gwwl) = "UCSC"
 g19 = liftOver( as(gwwl, "GRanges"), ch )
 message("liftover complete.")
 e = elementLengths(g19)
 dr = which(e==0)
 if (length(dr)>0) g19 = g19[-dr]
 g19 = unlist(g19)
 metadata(g19)$conversion = "liftOver"
 metadata(g19)$sessInfo.creation.liftOver = sessionInfo()
 seqlevelsStyle(g19) = "NCBI"
 genome(g19) = "GRCh37"
 cursn = seqlevels(g19)
 data(si.hs.37)
 seqinfo(g19) = si.hs.37[cursn]
 new("gwaswloc", extractDate=date(), g19)
}

