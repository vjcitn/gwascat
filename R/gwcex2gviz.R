#' 
#' Prepare salient components of GWAS catalog for rendering with Gviz
#' 
#' @importFrom S4Vectors mcols "mcols<-" metadata "metadata<-"
#' @importFrom GenomicFeatures exons
#' @importFrom AnnotationDbi mapIds
#' @importFrom IRanges overlapsAny
#' @param basegr gwaswloc instance containing information about GWAS in catalog
#' @param contextGR A GRanges instance delimiting the visualization in genomic
#' coordinates
#' @param txrefobj a TxDb instance
#' @param genesymobj an OrgDb instance
#' @param genome character tag like 'hg19'
#' @param plot.it logical, if FALSE, just return list
#' @param maxmlp maximum value of -10 log p -- winsorization of all larger
#' values is performed, modifying the contents of Pvalue\_mlogp in the
#' elementMetadata for the call
#' @keywords graphics
#' @examples
#' 
#'  data(ebicat_2020_04_30)
#'  # GenomeInfoDb::seqlevelsStyle(ebicat_2020_04_30) = "UCSC" # no more
#'  Seqinfo::seqlevels(ebicat_2020_04_30) = paste0("chr", Seqinfo::seqlevels(ebicat_2020_04_30))
#'  gwcex2gviz(ebicat_2020_04_30)
#' 
#' @export gwcex2gviz
gwcex2gviz = function( basegr, 
   contextGR = GRanges(seqnames="chr17", IRanges::IRanges(start=37500000, width=1e6)), 
   txrefobj = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, genome="hg19",
   genesymobj = org.Hs.eg.db::org.Hs.eg.db,
   plot.it=TRUE, maxmlp=25 ) {
#
# objective is to visualize features of GWAS in gene/transcript context
#
 if (!requireNamespace("Gviz")) stop("install Gviz to use this function")
 chrmin = as.character(seqnames(contextGR))
#
# the get() here is a hack.  need to have a way of getting relevant object
# name from txrefobj, indeed allowing more general spec of an exon range resource
#
 ex = exons( txrefobj, columns = c("gene_id", "tx_id", "exon_id"),
    filter=list(exon_chrom = chrmin) )
 txin = ex[ which(overlapsAny(ex, contextGR)) ]
 if (length(txin) == 0) stop("no transcripts in contextGR")
 v = mcols(txin)
 e = v$exon_id
 txl = v$tx_id
 texx = sapply(as.list(v$tx_id), "[", 1)
 g = sapply(as.list(v$gene_id), "[", 1)
 g = unlist(g)
 drop = which(is.na(g))
 k = GRanges(seqnames=chrmin, ranges=IRanges::ranges(txin), gene=g, exon=e,
    transcript=texx, id=1:length(g))
 if (length(drop) > 0) k = k[-drop]
 kk = mapIds(genesymobj, keys=mcols(k)$gene, keytype="ENTREZID",
           column="SYMBOL")  # multiVals == first
 mcols(k)$symbol = kk
 GR = Gviz::GeneRegionTrack(k, chromosome=chrmin, genome=genome)
 studs = basegr[ which(overlapsAny(basegr, contextGR)) ]
 mlp = mcols(studs)$PVALUE_MLOG
 mlp = ifelse(mlp > maxmlp, maxmlp, mlp)
 mcols(studs)$PVALUE_MLOG = mlp
 sss = mcols(studs)$PVALUE_MLOG
 studp = Gviz::DataTrack( as(studs, "GRanges"), data=sss, 
     chromosome=chrmin, genome=genome, name="-log P GWAS" )
#
# need to move these controls up to interface
#
 Gviz::displayPars(GR)$collapseTranscripts = TRUE
 Gviz::displayPars(GR)$showId = TRUE
 GR@name = "Genes"
 if (plot.it) Gviz::plotTracks(list(  studp, Gviz::GenomeAxisTrack(),  GR))
 invisible( list( studp, Gviz::GenomeAxisTrack(), GR ) )
}

