
intv = function(low,hi) {
  if (missing(low)) return("")
  paste0("[", low, "+TO+", hi, "]")
}
zpad = function(x){ cx = as.character(x); ifelse(nchar(cx)==1,
    paste0("0", cx), cx) }
datev = function(y0, m0, y1, m1) {
   stopifnot(nchar(y0)==4)
   paste0("[", y0, "-", zpad(m0), "-01T00:00:00Z+TO+", y1, "-",
      zpad(m1), "-01T00:00:00Z]")
}
tpad = function(tphr) gsub(" ", "+", phr)
tfilt = function(tphr) paste("&traitfilter[]=", tpad(phr), collapse="")

buildq = function(term, pvbound, orlim, betalim, y0, m0, y1, m1, 
  tphrvec, submit=TRUE ) {
if (missing(term)) stop("a search term must be provided")
if (missing(pvbound)) pvbound=""
if (missing(orlim)) orlimi=""
  else orlimi = intv(orlim[1], orlim[2])
if (missing(betalim)) betalimi=""
  else betalimi = intv(betalim[1], betalim[2])
if (missing(y0)) datei=""
  else datei = datev(y0,m0,y1,m1)
f1 = "http://www.ebi.ac.uk/gwas/api/search/downloads?q=text:"
f2 = "&pvalfilter="
f3 = "&orfilter="
f4 = "&betafilter="
f5 = "&datefilter="
f6 = "&dateaddedfilter="
oo = getOption("useFancyQuotes")

options(useFancyQuotes=FALSE)
on.exit(options(useFancyQuotes=oo))
equo = function(x) dQuote(x)

quer = paste0(f1, equo(term), f2, pvbound, f3, orlimi,
    f4, betalimi, f5, datei, f6)
if (!submit) return(quer)

suppressWarnings({
 read.delim(quer, sep="\t", header=TRUE, stringsAsFactors=FALSE)
})
}
