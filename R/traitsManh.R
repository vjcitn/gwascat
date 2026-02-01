#' use ggbio facilities to display GWAS results for selected traits in genomic
#' coordinates
#' 
#' use ggbio facilities to display GWAS results for selected traits in genomic
#' coordinates
#' 
#' uses a ggbio autoplot
#' 
#' @param gwr GRanges instance as managed by the gwaswloc container design,
#' with Disease.Trait and Pvalue\_mlog among elementMetadata columns
#' @param selr A GRanges instance to restrict the \code{gwr} for visualization.
#' Not tested for noncontiguous regions.
#' @param traits Character vector of traits to be exhibited; GWAS results with
#' traits not among these will be labeled ``other''.
#' @param truncmlp Maximum value of -log10 p to be displayed; in the raw data
#' this ranges to the hundreds and can cause bad compression.
#' @param \dots not currently used
#' @return autoplot value
#' @note An xlab is added, concatenating genome tag with seqnames tag.
#' @author VJ Carey <stvjc@@channing.harvard.edu>
#' @keywords models graphics
#' @examples
#' 
#' # do a p-value truncation if you want to reduce compression
#' \dontrun{  # ggbio July 2018
#' data(ebicat_2020_04_30)
#' library(GenomeInfoDb)
#' seqlevelsStyle(ebicat_2020_04_30) = "UCSC"
#' traitsManh(ebicat_2020_04_30)
#'  }
#' 
#' @export traitsManh
traitsManh = function( gwr, 
   selr=GRanges(seqnames="chr17", IRanges(3e7, 5e7)), 
   traits = c("Asthma", "Parkinson's disease", "Height", "Crohn's disease"),
   truncmlp = 25,
   ...)
 {
 Trait <- NA # try to squelch note
 if (!requireNamespace("ggbio")) stop("install ggbio to use this function")
 gwr = gwr[ which(overlapsAny(gwr, selr)) ]
 availtr = as.character(mcols(gwr)[["DISEASE/TRAIT"]])
 oth = which(!(availtr %in% traits))
 availtr[oth] = "Other"
 mcols(gwr)$Trait = availtr
 pv = mcols(gwr)$PVALUE_MLOG 
 mcols(gwr)$PVALUE_MLOG = ifelse(pv > 25, 25, pv)
 sn = paste(genome(gwr)[1], as.character(seqnames(gwr))[1], sep=" ")
 ggbio::autoplot(gwr, geom = "point", aes(y = PVALUE_MLOG, color = Trait),
      xlab=sn)

}
