gwcex2gviz = function( 
   contextGR = GRanges(seqnames="chr17", IRanges(start=37500000, width=1e6)), 
   txrefpk = "TxDb.Hsapiens.UCSC.hg19.knownGene", genome="hg19",
   genesympk = "org.Hs.eg.db",
   plot.it=TRUE ) {
#
# objective is to visualize features of GWAS in gene/transcript context
#
 require(Gviz, quietly=TRUE)
 library(txrefpk, character.only=TRUE, quietly=TRUE)
 library(genesympk, character.only=TRUE, quietly=TRUE)
 symmap = get(gsub(".db", "SYMBOL", genesympk))
 chrmin = as.character(seqnames(contextGR))
#
# the get() here is a hack.  need to have a way of getting relevant object
# name from txrefpk, indeed allowing more general spec of an exon range resource
#
 ex = exons( get(txrefpk), columns = c("gene_id", "tx_id", "exon_id"),
    vals=list(exon_chrom = chrmin) )
 txin = ex[ which(ex %in% contextGR ) ]
 if (length(txin) == 0) stop("no transcripts in contextGR")
 v = values(txin)
 e = v$exon_id
 txl = v$tx_id
 texx = sapply(as.list(v$tx_id), "[", 1)
 g = sapply(as.list(v$gene_id), "[", 1)
 g = unlist(g)
 drop = which(is.na(g))
 k = GRanges(seqnames=chrmin, ranges=ranges(txin), gene=g, exon=e,
    transcript=texx, id=1:length(g))
 if (length(drop) > 0) k = k[-drop]
 kk = unlist(mget(values(k)$gene, symmap))
 values(k)$symbol = kk
 GR = GeneRegionTrack(k, chromosome=chrmin, genome=genome)
 library(gwascat)
 studs = gwrngs[ which(gwrngs %in% contextGR) ]
 sss = values(studs)$Pvalue_mlog
 studp = DataTrack( as(studs, "GRanges"), data=sss, 
     chromosome=chrmin, genome=genome, name="-log P GWAS" )
#
# need to move these controls up to interface
#
 displayPars(GR)$collapseTranscripts = TRUE
 displayPars(GR)$showId = TRUE
 GR@name = "Genes"
 if (plot.it) plotTracks(list(  studp, GenomeAxisTrack(),  GR))
 invisible( list( studp, GenomeAxisTrack(), GR ) )
}

