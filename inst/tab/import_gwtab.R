import_gwtab = function(fn) {
 #
 # shouldn't be necessary but there are odd phenomena in the
 # distributed file
 #
 recs = readLines(fn)
 nc = nchar(recs)
 if (any(nc==0)) recs = recs[-which(nc==0)] # shld drop last
 fhead = strsplit(recs[1], "\t")[[1]]
 dat = strsplit(recs[-1], "\t")
 nrec = length(dat)
 nfields = length(fhead)
 nfdat = sapply(dat, length)
 if (length(table(nfdat)) > 1) stop("irregular field number in records, sorry.")
 outdat = matrix("", nr=nrec, nc=nfields)
 for (i in 1:nrec) outdat[i,] = dat[[i]]
 colnames(outdat) = fhead
 data.frame(outdat, stringsAsFactors=FALSE, check.names=FALSE)
}
