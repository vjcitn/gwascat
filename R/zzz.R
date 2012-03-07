.onAttach = function(libname, pkgname) {
#
# create global data objects
#  1) gwcat, a data.frame instance directly reflecting content of the table from NHGRI
#  2) gwrngs, a GRanges that is filtered to studies with specific claims of SNP-trait associations
#
 gwcat <<- get(load(system.file("data/gwdf_2012_03_07.rda", package="gwascat")))
#
# !! please reset extractDate as appropriate
#
 extractDate = "2012.03.07"
 psm =  function(..., appendLF=FALSE )packageStartupMessage(..., appendLF=appendLF)
 psm(paste("'gwcat' data frame now available, provides NHGRI GWAS cat records of ", extractDate,".\n", sep=""))
 psm("building 'gwrngs', GRanges for studies with located variants...")
 gwcatloc = gwcat[nchar(gwcat$Chr_pos)>0,]
#
# various fixes
# 1) put chrnn
#
 ch = gwcatloc$Chr_id
 if (length(grep("chr", ch)) == 0) ch = paste("chr", ch, sep="")
 gwrngs = GRanges(seqnames=ch, IRanges(as.numeric(gwcatloc$Chr_pos),width=1))
 values(gwrngs) = gwcatloc
#
# 2) make numeric p values and addresses
#
 values(gwrngs)$p.Value = as.numeric(values(gwrngs)$p.Value)
 values(gwrngs)$Chr_pos = as.numeric(values(gwrngs)$Chr_pos)
#
# 3) clean out stray whitespace
#
 badco = values(gwrngs)$"Strongest.SNP.Risk.Allele"
 co = gsub(" $", "", badco)
 values(gwrngs)$"Strongest.SNP.Risk.Allele" = co
#
# 4) deal with OR or beta field entries possessing strings
#
 fixhet = function(vec) {
# take a mix of numerical strings and character strings 
# and set char strings to ""
# so that as.numeric will not warn
   strinds = grep("[a-zA-Z]", vec)
   if (length(strinds)>0) vec[strinds] = ""
   vec
   }
# strinds = grep("[a-zA-Z]", values(gwrngs)[,"OR.or.beta"])
# if (length(strinds)>0) values(gwrngs)$OR.or.beta[strinds] = ""
 values(gwrngs)$OR.or.beta = as.numeric(fixhet(values(gwrngs)$OR.or.beta))
 values(gwrngs)$Risk.Allele.Frequency = as.numeric(fixhet(values(gwrngs)$Risk.Allele.Frequency))
 gwrngs = new("gwaswloc", extractDate=extractDate, gwrngs)
 assign("gwrngs", gwrngs, .GlobalEnv)
 psm("done.\n")
}
