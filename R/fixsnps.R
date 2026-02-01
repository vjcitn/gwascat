#' this function looks for clues to resolve records that have missing CHR_ID or CHR_POS,
#' for which address info can be inferred from SNPS field
#' @param x output of makeCurrentGwascat or get_cached_gwascat
#' @export
fixsnps = function(x) {
 nochrid = which(nchar(x[["CHR_ID"]])==0)
 nochrpos = which(nchar(x[["CHR_POS"]])==0)
 badrecs = unique(c(nochrid, nochrpos))
 if (length(badrecs)==0) return(x)
 cands = x[badrecs, "SNPS"][[1]] # convert to char
 isrs = grep("^rs", cands)
 if (length(isrs)>0) {
   cands = cands[-isrs]
   badrecs = badrecs[-isrs]
   }
 scands = strsplit(cands, ":")
 idcands = sapply(scands, "[", 1)
 poscands = sapply(scands, "[", 2)
 idcands = gsub("chr", "", idcands)
 ok = which(idcands %in% as.character(c(1:22, "X", "Y", "MT")))
 badrecs = badrecs[ok]
 poscands = poscands[ok]
 idcands = idcands[ok]
 notnum = suppressWarnings(which(is.na(as.numeric(poscands))))
 if (length(notnum)>0) {
  badrecs = badrecs[-notnum]
  poscands = poscands[-notnum]
  idcands = idcands[-notnum]
  }
 x[badrecs, "CHR_ID"] = idcands
 x[badrecs, "CHR_POS"] = poscands
 ok = which(x$CHR_ID %in% as.character(c(1:22, "X", "Y", "MT")) & suppressWarnings(as.numeric(x$CHR_POS))<(2^31-1))
 noncanon = x[-ok,]
 x = x[ok,]
 attr(x, "noncanon") = noncanon
 x
}
