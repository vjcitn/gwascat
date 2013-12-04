
makeCurrentGwascat = function(table.url="http://www.genome.gov/admin/gwascatalog.txt", fixNonASCII=TRUE) {
 tab = read.delim(url(table.url), sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
 if (fixNonASCII) tab = fixNonASCII(tab)
 gwdf2GRanges(tab, extractDate=as.character(Sys.Date()))
}
