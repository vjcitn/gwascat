traitsManh = function( gwr, 
   selr=GRanges(seqnames="chr17", IRanges(3e7, 5e7)), 
   traits = c("Asthma", "Parkinson's disease", "Height", "Crohn's disease"),
   truncmlp = 25,
   ...)
 {
 Trait <- NA # try to squelch note
 require(ggbio)
 gwr = gwr[ which(overlapsAny(gwr, selr)) ]
 availtr = as.character(mcols(gwr)$DISEASE.TRAIT)
 oth = which(!(availtr %in% traits))
 availtr[oth] = "Other"
 mcols(gwr)$Trait = availtr
 pv = mcols(gwr)$PVALUE_MLOG 
 mcols(gwr)$PVALUE_MLOG = ifelse(pv > 25, 25, pv)
 autoplot( gwr, geom="point", aes(y=PVALUE_MLOG, color=Trait))
}
