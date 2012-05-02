
makeCurrentGwascat = function(table.url="http://www.genome.gov/admin/gwascatalog.txt", fixNonASCII=TRUE) {
 tab = read.delim(url(table.url), sep="\t", header=TRUE)
 if (fixNonASCII) tab = gwascat:::fixNonASCII(tab)
 gwascat:::gwdf2GRanges(tab, extractDate=as.character(Sys.Date()))
}
