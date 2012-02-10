.onAttach = function(libname, pkgname) {
#
# create global data objects
#  1) gwcat, a data.frame instance directly reflecting content of the table from NHGRI
#  2) gwrngs, a GRanges that is filtered to studies with specific claims of SNP-trait associations
#
 gwcat <<- get(load(system.file("data/gwdf_2012_02_02.rda", package="gwascat")))
 psm =  function(..., appendLF=FALSE )packageStartupMessage(..., appendLF=appendLF)
 psm("'gwcat' data frame now available, provides NHGRI GWAS cat records of 02/02/2012.\n")
 psm("building 'gwrngs', GRanges for studies with located variants...")
 gwcatloc = gwcat[nchar(gwcat$Chr_pos)>0,]
 ch = gwcatloc$Chr_id
 if (length(grep("chr", ch)) == 0) ch = paste("chr", ch, sep="")
 gwrngs = GRanges(seqnames=ch, IRanges(as.numeric(gwcatloc$Chr_pos),width=1))
 values(gwrngs) = gwcatloc
 values(gwrngs)$p.Value = as.numeric(values(gwrngs)$p.Value)
 values(gwrngs)$Chr_pos = as.numeric(values(gwrngs)$Chr_pos)
 strinds = grep("[A-Z]", values(gwrngs)[,"OR.or.beta"])
 if (length(strinds)>0) values(gwrngs)$OR.or.beta[strinds] = ""
 values(gwrngs)$OR.or.beta = as.numeric(values(gwrngs)$OR.or.beta)
 gwrngs = new("gwaswloc", gwrngs)
 assign("gwrngs", gwrngs, .GlobalEnv)
 psm("done.\n")
}
