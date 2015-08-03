
topTraits = function(gwwl, n=10, tag="DISEASE.TRAIT") {
 sort(table(mcols(gwwl)[[tag]]), decreasing=TRUE)[1:n]
}

locs4trait = function(gwwl, trait, tag="DISEASE.TRAIT") {
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
  require(SNPlocs.Hsapiens.dbSNP.20120608)
  allrs = mcols(gwwl)$SNPS
  allch = mcols(gwwl)$CHR_ID
  rsbyc = split(allrs, allch)
  rcur = rsbyc[[chrtag]]
  refcur = getSNPlocs(paste("ch", chrtag, sep=""))
  rownames(refcur) = paste("rs", refcur$RefSNP_id, sep="")
  Ncur = refcur[ intersect(rcur, rownames(refcur)), ]
  minds <- match(rownames(Ncur), mcols(gwwl)$SNPS)
  #Ncurpos = mcols(gwwl[ minds ])$CHR_POS
  Ncurpos = start(gwwl[ minds ])
  all((Ncurpos - Ncur[,3]) == 0)
}

variantProps = function(rs, ..., gwwl=gwrngs) {
  subr = gwwl[rs,]
  strg = mcols(subr)$STRONGEST.SNP.RISK.ALLELE
  alls = sapply(strsplit(strg, "-"), "[", 2)
  strs = sapply(strsplit(strg, "-"), "[", 1)
  mcols(subr) = DataFrame(rsid=strs, riskAllele=alls,
     mcols(subr)[, c("DISEASE.TRAIT", "SNPS", "P.VALUE")])
  subr
}



snpGenos = function(chr, snpap="SNPlocs.Hsapiens.dbSNP.20120608") {
#
# for S snp on chrom chr, returns an S x 3 matrix with columns A/A, A/B, B/B
# the i, j element is the nuc assignment for genotype j of snp i
#
 require(snpap, character.only=TRUE, quietly=TRUE)
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

AB2nuc = function(abvec, rsids, chr, snpannopk="SNPlocs.Hsapiens.dbSNP.20120608") {
#
# converts variant calls in A/B form to nucleotides assuming A in A/B is
# alphabetically earlier and nucleotide mappings defined for rsids are as in
# snpannopk
#
 require(snpannopk, character.only=TRUE, quietly=TRUE)
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

 
ABmat2nuc = function(abmat, chr, snpannopk="SNPlocs.Hsapiens.dbSNP.20120608",
  gencode = c("A/A", "A/B", "B/B")) {
#
# converts variant calls in A/B form to nucleotides assuming A in A/B is
# alphabetically earlier and nucleotide mappings defined for rsids are as in
# snpannopk
#
 require(snpannopk, character.only=TRUE, quietly=TRUE)
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
   gwwl, snpap="SNPlocs.Hsapiens.dbSNP.20120608",
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
 
 
 
