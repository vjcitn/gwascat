
addSeqlengths = function( gr, src = "TxDb.Hsapiens.UCSC.hg19.knownGene" ) {
 require(src, character.only=TRUE)
 tx = transcripts(get(src))
 sl = seqlengths(tx)
 sl = sl[ names(seqlengths(gr)) ]
 seqlengths(gr) = sl
 gr
}
