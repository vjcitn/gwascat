
topTraits = function(gwwl, n=10, tag="DISEASE/TRAIT") {
 sort(table(mcols(gwwl)[[tag]]), decreasing=TRUE)[1:n]
}

locs4trait = function(gwwl, trait, tag="DISEASE/TRAIT") {
 if (length(trait) != 1) stop("please supply only one trait string")
 curem = mcols(gwwl) #as(mcols(gwwl), "DataFrame")
 tr = curem[[tag]]
 if (!(trait %in% tr)) stop(paste(trait, "not found in ", substitute(gwwl)))
 okind = which(tr == trait)
 new("gwaswloc", gwwl[okind])
}

chklocs = function(chrtag="20", gwwl=gwrngs19) {
#
# return TRUE if all named SNPs with locations in both
# the SNPlocs package and the gwascat agree 
#
# Note -- as of 2014, we cannot use the Chr_pos element of
# the shipped catalog.  We have to use the result of liftOver,
# in the gwrngs19 structure
#
  requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37")
  locdata = SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37
  allrs = mcols(gwwl)$SNPS
  allch = mcols(gwwl)$CHR_ID # is NCBI format, we hope
  rsbyc = split(allrs, allch)
  rcur = rsbyc[[chrtag]]
  refcur = snpsBySeqname(locdata, chrtag) #GPos with rsid in $RefSNP_id
  ref2keep_inds = match(rcur, refcur$RefSNP_id, nomatch=0)
  chkstat = setdiff(start(gwwl[rcur[which(ref2keep_inds != 0)]]) ,
            start(refcur[ref2keep_inds]))
  if (length(chkstat) > 0) {
       warning(paste0(length(chkstat), " snp ids in gwascat had reported location discordant with reference address; printing:"), immediate.=TRUE)
       print(gwwl[ which(start(gwwl[rcur[which(ref2keep_inds != 0)]]) %in% chkstat) ])
       return(FALSE)
  }
  return(TRUE)
}

variantProps = function(rs, ..., gwwl=gwrngs) {
  subr = gwwl[rs,]
  strg = mcols(subr)[["STRONGEST SNP-RISK ALLELE"]]
  alls = sapply(strsplit(strg, "-"), "[", 2)
  strs = sapply(strsplit(strg, "-"), "[", 1)
  mcols(subr) = DataFrame(rsid=strs, riskAllele=alls,
     mcols(subr)[, c("DISEASE/TRAIT", "SNPS", "P-VALUE")])
  subr
}



snpGenos = function(chr, snpap="SNPlocs.Hsapiens.dbSNP144.GRCh37") {
#
# for S snp on chrom chr, returns an S x 3 matrix with columns A/A, A/B, B/B
# the i, j element is the nuc assignment for genotype j of snp i
#
 requireNamespace(snpap, quietly=TRUE)
 df = getSNPlocs(chr)
 rsids = paste("rs", df$RefSNP_id, sep="")
 alls = as.character(df$alleles_as_ambig)
 map = Biostrings::IUPAC_CODE_MAP
 lets = map[alls]
 hom1 = sub("(.)(.)", "\\1/\\1", lets)
 hom2 = sub("(.)(.)", "\\2/\\2", lets)
 het = sub("(.)(.)", "\\1/\\2", lets)
 ans = cbind(hom1, het, hom2)
 rownames(ans) = rsids
 colnames(ans) = c("A/A", "A/B", "B/B")
 ans
}

AB2nuc = function(abvec, rsids, chr, snpannopk="SNPlocs.Hsapiens.dbSNP144.GRCh37") {
#
# converts variant calls in A/B form to nucleotides assuming A in A/B is
# alphabetically earlier and nucleotide mappings defined for rsids are as in
# snpannopk
#
 requireNamespace(snpannopk, quietly=TRUE)
 if (length(chr) > 1) stop("supports only chr of length 1")
#
# process Bioc SNP location metadata
#
 mapmat = snpGenos( chr, snpap = snpannopk )
 refrsid = rownames(mapmat)
#
# process input codes and ids
#
 names(abvec) = rsids
 feasible_ids = intersect(refrsid, rsids)
 nlost = length(setdiff(rsids, refrsid))
 abvec = abvec[feasible_ids]
 mapmat = mapmat[feasible_ids,]
 bad = is.na(abvec) | !(abvec %in% c("A/A", "A/B", "B/B"))
 if (any(bad)) {
    abvec = abvec[-which(bad)]
    feasible_ids = feasible_ids[-which(bad)]
    mapmat = mapmat[feasible_ids,]
    }
#
# now just pull the columns corresponding to the A/B found
#
 ans = mapmat[ cbind(feasible_ids, abvec) ]
 names(ans) = feasible_ids
 ans
}

 
ABmat2nuc = function(abmat, chr, snpannopk="SNPlocs.Hsapiens.dbSNP144.GRCh37",
  gencode = c("A/A", "A/B", "B/B")) {
#
# converts variant calls in A/B form to nucleotides assuming A in A/B is
# alphabetically earlier and nucleotide mappings defined for rsids are as in
# snpannopk
#
 requireNamespace(snpannopk, quietly=TRUE)
 if (length(chr) > 1) stop("supports only chr of length 1")
#
# process Bioc SNP location metadata
#
 mapmat = snpGenos( chr, snpap = snpannopk )
 refrsid = rownames(mapmat)
#
# process input codes and ids
#
 rsids = colnames(abmat)
 feasible_ids = intersect(refrsid, rsids)
 nlost = length(setdiff(rsids, refrsid))
 abmat = abmat[,feasible_ids]
 mapmat = mapmat[feasible_ids,]
 if (!all(abmat %in% gencode)) {
   abmat[ !(abmat %in% gencode) ] = "N/N"
   mapmat = cbind(mapmat, "N/N")
   colnames(mapmat)[4] = "N/N"
 }
#
# now we have nrow(abmat) rows to process
#
#
# now just pull the columns corresponding to the A/B found
#
 ans = apply(abmat, 1, function(x) 
        mapmat[ cbind(feasible_ids, x) ] )
 rownames(ans) = feasible_ids
 t(ans)
}

riskyAlleleCount = function(callmat, matIsAB=TRUE, chr,
   gwwl, snpap="SNPlocs.Hsapiens.dbSNP144.GRCh37",
   gencode = c("A/A", "A/B", "B/B")) {
 if (matIsAB) callmat = ABmat2nuc( abmat=callmat, chr=chr, snpannopk=snpap, gencode=gencode)
 uchr = gsub("ch", "chr", chr)
 gwwl = subsetByChromosome(gwwl, uchr)
 possrs = getRsids(gwwl)
 callmat = callmat[, intersect(colnames(callmat), possrs)]
 possrs = colnames(callmat)
 gwwl = gwwl[possrs,]
 vp = variantProps(possrs, gwwl=gwwl)
 list(callmat=callmat, vp=vp)
 vpallele = mcols(vp)$riskAllele
 vpallele = gsub("\\?", "@@", vpallele)
 mcols(vp)$riskAllele = vpallele
 nhits = t(apply(callmat, 1, function(x) {
      sapply(1:length(x), function(z) length(grep(vpallele[z],
               strsplit(x[z], "/")[[1]]))) }))
 colnames(nhits) = possrs
 nhits
}
 
 
 
uri2node = function(us) {
 stopifnot(length(us)==1, is.atomic(us))
 ents = strsplit(us, ", ")[[1]]
 gsub("http://www.ebi.ac.uk/efo/EFO_", "EFO:", ents)
 }

node2uri = function(nn) {
 stopifnot(length(nn)==1, is.atomic(nn))
 paste0("http://www.ebi.ac.uk/efo/", gsub(":", "_",nn))
 }
