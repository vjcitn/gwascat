#' convert GWAS catalog data.frame to gwaswloc, a GRanges extension with simple show method
#' @param df data.frame
#' @export
process_gwas_dataframe = function (df) {
#
# intent is to take a data frame like that distributed by EMBL/EBI (formerly by NHGRI)
# and convert to a useful GRanges instance, coercing heterogeneous vectors
# to majority type
#
# NOTE: EMBL/EBI changed the column header case to capitals.  Code
#  will use the new naming convention
#
# April 30 2020 -- we find that CHR_POS can include NA and c x c and c X c for SNP-SNP interaction studies
# We cannot make these into GRanges.  We hand them back in metadata component
    badpos = unique(c(
              which(is.na(df$CHR_ID)), 
              which(is.na(df$CHR_POS)), grep(";|x|X", df$CHR_POS), which(df$CHR_POS=="")))
    nogr = df[badpos,]
    gwcatloc = df[-badpos,]
#
# put chr prefix to CHR_ID as needed
#
    ch = as.character(gwcatloc$CHR_ID)
#    if (any(ch == "23")) ch[ch=="23"] = "X"
    
    if (length(grep("chr", ch)) == 0) 
        ch = paste("chr", ch, sep = "")
    gwrngs = GRanges(seqnames = ch, IRanges(as.numeric(gwcatloc$CHR_POS), 
        width = 1))
    mcols(gwrngs) = gwcatloc
#
# make numeric p values and addresses
#
    mcols(gwrngs)[["P-VALUE"]] = 
  as.numeric(as.character(mcols(gwrngs)[["P-VALUE"]])) # was factor
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
# utility to get numeric values in Risk.Allele.Frequency -- note, there may be scientific notation
#
    killpatt = "\\-|\\+|[[:alpha:]]|\\(|\\)|\\ "
#    nulcToNA = function(x) {isn = which(nchar(x)==0); if (length(isn)>0) x[isn] = NA; x}
    raf = mcols(gwrngs)[["RISK ALLELE FREQUENCY"]] 
    suppressWarnings({raf = as.numeric(raf)}) # will make NA for 'NR', ranges, etc. but keep scientific notation
    mcols(gwrngs)[["RISK ALLELE FREQUENCY"]] = raf
    gwrngs = new("gwaswloc", extractDate = attr(df, "extractDate"), gwrngs)
    GenomeInfoDb::seqlevelsStyle(gwrngs) = "NCBI"
    data("si.hs.38", package="gwascat")
    sl = seqlevels(gwrngs)
    Seqinfo::seqinfo(gwrngs) = si.hs.38[sl]
    metadata(gwrngs) = list(
      date.created = date(),
      creation = match.call(),
      badpos = nogr,
      sessInfo.creation = sessionInfo()
    )
    gwrngs
}
