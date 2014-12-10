
ldtagr = function( snprng, tf, samples, genome="hg19",
   lbmaf=.05, lbR2=.8, radius=100000 ) {
#
# given a GRanges with name and address of a SNP with representation in a
# VCF referenced by tf, give the names and addresses of all
# SNP in LD with the original locus at R2 at least lbR2
# 
  stopifnot(length(snprng)==1)
  snpid = names(snprng)
  stopifnot(length(snpid)==1)
  quer = gQTLstats:::queryVCF( gr=snprng+radius, vcf.tf=tf, 
         samps=samples, genome=genome, getSM=TRUE )
  empty = GRanges()
  mcols(empty) = DataFrame(paramRangeID = factor(), R2 = numeric())
  vcfrng = rowData(quer$readout)
  gt = quer$sm$genotypes
  if (!(snpid %in% colnames(gt))) {
    message(paste0("NOTE: ", snpid, " not in VCF at given radius, returning empty GRanges"))
    return(empty)
    }
  map = quer$sm$map
  cs = col.summary(gt)
  mafok = which(cs[,"MAF"] >= lbmaf)
  stopifnot(length(mafok)>0)
  gt = gt[,mafok]
  if (!(snpid %in% colnames(gt))) {
    message(paste0("NOTE: ", snpid, " has MAF too low, returning empty GRanges."))
    return(empty)
    }
  ldvec = ld(gt[,snpid], gt, stats="R.squared")
  extrng = vcfrng[ thec <- colnames(theld <- ldvec[, which(as.numeric(ldvec) > lbR2),drop=FALSE ]) ]
  extrng$R2 = theld
  #list(gt=gt, map=map, vcfrng=vcfrng, ldvec=ldvec,
  #   extrng = vcfrng[ colnames(ldvec[, 
  #   which(as.numeric(ldvec) > lbR2),drop=FALSE ]) ])
  names(extrng) = make.names(thec)
  extrng[, c("paramRangeID", "R2")]
}
