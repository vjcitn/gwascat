#' use AnnotationHub snapshot as basis for gwaswloc structure creation
#' @param x inherits from data.frame, with columns consistent with EBI table
#' @param fixNonASCII logical(1) if TRUE, use iconv to replace non-ASCII character, important
#' for CMD check but perhaps not important for applied use
#' @examples
#' ah = AnnotationHub::AnnotationHub()
#' entitytab = AnnotationHub::query(ah, "gwascatData")
#' cand = names(entitytab)[1]
#' stopifnot(nchar(cand)>0)
#' tab = ah[[cand]]
#' gww = gwcat_snapshot(tab)
#' gww
#' length(gww)
#' @export
gwcat_snapshot = function(x, fixNonASCII=TRUE) {
#    ct = readr::cols(.default = col_character(), `DATE ADDED TO CATALOG` = col_date(format = ""), 
#        PUBMEDID = col_double(), DATE = col_date(format = ""), 
#        CHR_ID = col_character(), CHR_POS = col_double(), UPSTREAM_GENE_DISTANCE = col_double(), 
#        DOWNSTREAM_GENE_DISTANCE = col_double(), MERGED = col_double(), 
#        SNP_ID_CURRENT = col_double(), INTERGENIC = col_double(), 
#        `P-VALUE` = col_double(), PVALUE_MLOG = col_double(), 
#        `OR or BETA` = col_double())
#    tab = readr::read_tsv(tf, col_types = ct)
#    message(paste0("formatting gwaswloc instance..."))
    tab = as.data.frame(x)
    if (fixNonASCII) 
        tab = fixNonASCII(tab)
    cur_plus = gwdf2GRanges(tab, extractDate = as.character(Sys.Date()))
    cur = cur_plus$okrngs
    nogr = cur_plus$nogr
    seqlevelsStyle(cur) = "NCBI"
    cursn = seqlevels(cur)
    data(si.hs.38)
    seqinfo(cur) = si.hs.38[cursn]
    metadata(cur) = list(date.created = "2021-03-30", creation = "AnnotationHub", 
        badpos = nogr, sessInfo.creation = sessionInfo())
    message("done.")
    cur
}

