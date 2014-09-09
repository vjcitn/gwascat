
addSeqlengths = function( gr, src ) {
 stopifnot(inherits(src, "Seqinfo"))
 sl = seqlengths(src)
 sl = sl[ names(seqlengths(gr)) ]  # seqlengths(gr) will be named even if NA
 seqlengths(gr) = sl    # should probably use seqinfo(gr)
 gr
}
