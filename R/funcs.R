
topTraits = function(gwwl, n=10, tag="Disease.Trait") {
 sort(table(values(gwwl)[[tag]]), decreasing=TRUE)[1:n]
}

locs4trait = function(gwwl, trait, tag="Disease.Trait") {
 if (length(trait) != 1) stop("please supply only one trait string")
 curem = elementMetadata(gwwl) #as(values(gwwl), "DataFrame")
 tr = curem[[tag]]
 if (!(trait %in% tr)) stop(paste(trait, "not found in ", substitute(gwwl)))
 okind = which(tr == trait)
 new("gwaswloc", gwwl[okind])
}

chklocs = function(chrtag="20", gwwl=gwrngs) {
#
# return TRUE if all named SNPs with locations in both
# the SNPlocs package and the gwascat agree 
#
  require(SNPlocs.Hsapiens.dbSNP.20110815)
  allrs = elementMetadata(gwwl)$SNPs
  allch = elementMetadata(gwwl)$Chr_id
  rsbyc = split(allrs, allch)
  rcur = rsbyc[[chrtag]]
  refcur = getSNPlocs(paste("ch", chrtag, sep=""))
  rownames(refcur) = paste("rs", refcur$RefSNP_id, sep="")
  Ncur = refcur[ intersect(rcur, rownames(refcur)), ]
  minds <- match(rownames(Ncur), elementMetadata(gwwl)$SNPs)
  Ncurpos = elementMetadata(gwwl[ minds ])$Chr_pos
  all((as.numeric(elementMetadata(gwwl[ minds ])$Chr_pos) - Ncur[,3]) == 0)
}

