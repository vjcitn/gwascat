 fixhet = function(vec) {
# take a mix of numerical strings and character strings 
# and set char strings to ""
# so that as.numeric will not warn
   strinds = grep("[a-zA-Z]", vec)
   if (length(strinds)>0) vec[strinds] = ""
   vec
   }
.onAttach = function(libname, pkgname) {
packageStartupMessage("gwascat loaded.  Use data(ebicat38) for hg38 coordinates;")
packageStartupMessage(" data(ebicat37) for hg19 coordinates.")
}

#.onAttach = function(libname, pkgname) {
##
## create global data objects
##  1) gwcat, a data.frame instance directly reflecting content of the table from NHGRI
##  2) gwrngs, a GRanges that is filtered to studies with specific claims of SNP-trait associations
##
## gwcat <- get(load(system.file("data/gwdf_2012_09_22.rda", package="gwascat")))
## gwcat <<- fixNonASCII(gwcat)
##
## !! please reset extractDate as appropriate
##
## extractDate = "2013.12.03"
### psm =  function(..., appendLF=FALSE )packageStartupMessage(..., appendLF=appendLF)
## psm(paste("'gwcat' data frame now available, provides NHGRI GWAS cat records of ", extractDate,".\n", sep=""))
##if (0) {
## psm("building 'gwrngs', GRanges for studies with located variants...", appendLF=TRUE)
## gwcatloc = gwcat[nchar(gwcat$Chr_pos)>0,]
## assign("gwrngs", gwdf2GRanges(gwcat, extractDate=extractDate), .GlobalEnv)
##}
#psm =  function(..., appendLF=FALSE )packageStartupMessage(..., appendLF=appendLF)
#gwrngs <<- get(load(system.file("data/gwrngs.rda", package="gwascat")))
#psm("Object 'gwrngs' loaded and assigned from serialized version of 2013.12.03.", appendLF=TRUE)
#psm("Use makeCurrentGwascat() to obtain up-to-date image.", appendLF=TRUE)
#}

gwdf2GRanges = function (df, extractDate, seqlSrc ) 
{
#
# intent is to take a data frame like that distributed by EMBL/EBI (formerly by NHGRI)
# and convert to a useful GRanges instance, coercing heterogeneous vectors
# to majority type
#
# NOTE: EMBL/EBI changed the column header case to capitals.  Code
#  will use the new naming convention
#
    gwcatloc = df[which(!is.na(as.numeric(df$CHR_POS))), ]
#
# put chr prefix to CHR_ID as needed
#
    ch = as.character(gwcatloc$CHR_ID)
    if (any(ch == "23")) ch[ch=="23"] = "X"
    
    if (length(grep("chr", ch)) == 0) 
        ch = paste("chr", ch, sep = "")
    gwrngs = GRanges(seqnames = ch, IRanges(as.numeric(gwcatloc$CHR_POS), 
        width = 1))
    mcols(gwrngs) = gwcatloc
#
# make numeric p values and addresses
#
    mcols(gwrngs)$P-VALUE = as.numeric(as.character(mcols(gwrngs)$P.VALUE)) # was factor
    mcols(gwrngs)$PVALUE_MLOG = as.numeric(as.character(mcols(gwrngs)$PVALUE_MLOG)) # was factor
    mcols(gwrngs)[["OR or BETA"]] = suppressWarnings(as.numeric(as.character(mcols(gwrngs)[["OR or BETA"]]))) # was factor
    mcols(gwrngs)$CHR_POS = as.numeric(mcols(gwrngs)$CHR_POS)
#
# clean out stray whitespace
#
    badco = mcols(gwrngs)[["STRONGEST SNP-RISK ALLELE"]]
    co = gsub(" $", "", badco)
    mcols(gwrngs)[["STRONGEST SNP-RISK ALLELE"]] = co
#
#
#
#    fixhet = function(vec) {
#        strinds = grep("[a-zA-Z]", vec)
#        if (length(strinds) > 0) 
#            vec[strinds] = ""
#        vec
#    }
#
# utility to strip out some unusual tokens in OR.or.BETA
#
#    cleanToNum = function(x, bad = c("NR", "Pending")) {
#      lkbad = which(x %in% bad)
#      if (length(lkbad) > 0) x[lkbad] = NA
#    }
#    mcols(gwrngs)[["OR or BETA"]] = cleanToNum(
#          mcols(gwrngs)[["OR or BETA"]] ) #as.numeric(fixhet(mcols(gwrngs)[["OR or BETA"]]))
#
# utility to get numeric values in Risk.Allele.Frequency
#
    killpatt = "\\+|[[:alpha:]]|\\(|\\)|\\ "
    nulcToNA = function(x) {isn = which(nchar(x)==0); if (length(isn)>0) x[isn] = NA; x}
    mcols(gwrngs)[["RISK ALLELE FREQUENCY"]] = as.numeric(nulcToNA(gsub(killpatt, "", as.character(mcols(gwrngs)[["RISK ALLELE FREQUENCY"]]))))
#    gwrngs = makeConsecChrs(gwrngs)  # WHY???
#    gwrngs = addSeqlengths(gwrngs, src=seqlSrc)
    gwrngs = new("gwaswloc", extractDate = extractDate, gwrngs)
    gwrngs
}

