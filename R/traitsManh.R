traitsManh = function( gwr, 
   selr=GRanges(seqnames="chr17", IRanges(3e7, 5e7)), 
   traits = c("Asthma", "Parkinson's disease", "Height", "Crohn's disease"),
   truncmlp = 25,
   ...)
 {
 require(ggbio)
 gwr = gwr[ which(overlapsAny(gwr, selr)) ]
 availtr = mcols(gwr)$Disease.Trait
 oth = which(!(availtr %in% traits))
 availtr[oth] = "Other"
 mcols(gwr)$Trait = availtr
 pv = mcols(gwr)$Pvalue_mlog 
 mcols(gwr)$Pvalue_mlog = ifelse(pv > 25, 25, pv)
 autoplot( gwr, geom="point", aes(y=Pvalue_mlog, color=Trait))
}
