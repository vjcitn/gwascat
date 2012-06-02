
makeConsecChrs = function( gr ) {
 # assumes that chr names are lexically consecutive, want numerically consec.
 # e.g., chr1, chr10, chr11, ... should be chr1, chr2, chr3 ... mainly for use with ggbio::layout_circle
 fixvec = c(1,12,16,17,18,19,20,21,22,2,3,4,5,6,7,8,9,10,11,13,14,15,23)
 tmpinf = seqinfo( gr )
 tmpinf@seqnames = tmpinf@seqnames[fixvec]
 tmpinf@seqlengths = tmpinf@seqlengths[fixvec]
 seqlevels(gr) = seqlevels(gr)[fixvec]
 seqinfo(gr) = tmpinf
 gr
}

