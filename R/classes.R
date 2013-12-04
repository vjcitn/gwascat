
#> names(gwcat)
# [1] "Date Added to Catalog"      "PUBMEDID"                  
# [3] "First Author"               "Date"                      
# [5] "Journal"                    "Link"                      
# [7] "Study"                      "Disease/Trait"             
# [9] "Initial Sample Size"        "Replication Sample Size"   
#[11] "Region"                     "Chr_id"                    
#[13] "Chr_pos"                    "Reported Gene(s)"          
#[15] "Mapped_gene"                "Upstream_gene_id"          
#[17] "Downstream_gene_id"         "Snp_gene_ids"              
#[19] "Upstream_gene_distance"     "Downstream_gene_distance"  
#[21] "Strongest SNP-Risk Allele"  "SNPs"                      
#[23] "Merged"                     "Snp_id_current"            
#[25] "Context"                    "Intergenic"                
#[27] "Risk Allele Frequency"      "p-Value"                   
#[29] "Pvalue_mlog"                "p-Value (text)"            
#[31] "OR or beta"                 "95% CI (text)"             
#[33] "Platform [SNPs passing QC]" "CNV"              
#
#setClass("study", representation(author="character",
# pmid="character", date="character"

setClass("gwaswloc", representation(extractDate="character"), 
   contains="GRanges")

setMethod("show", "gwaswloc", function(object) {
 cat("gwasloc instance with", length(object), "records and", ncol(mcols(object)),
  "attributes per record.\n")
 cat("Extracted: ", object@extractDate, "\n")
 cat("Excerpt:\n")
 nrec = min(5, length(object))
 availcols = colnames(mcols(object))
 majorfields = c("Disease.Trait", "SNPs", "p.Value")
 if (all(majorfields %in% availcols) & length(availcols) > 5)
   show(as(object, "GRanges")[1:nrec, majorfields])
 else show(as(object, "GRanges")[1:nrec,])
})

#
# intention here is to allow row subscripting by rs number without any
# casting ... and to allow all other methods to proceed as in GRanges API
#
setMethod("[", "gwaswloc", function(x, i, j, ..., drop=FALSE) {
 if (missing(drop)) drop <- FALSE
# if (!missing(j)) stop("no column subscripting on gwaswloc")
 if (missing(i)) i = 1:length(x)
 if (!is(i, "numeric")) {
  rsids = getRsids(x)
  i = match(i, rsids, nomatch=0)
 }
 if (length(i)==0) stop("index has length 0")
 callNextMethod()
})
 

setGeneric("getRsids", function(x)standardGeneric("getRsids"))
setMethod("getRsids", "gwaswloc", function(x)
 mcols(x)$SNPs)

setGeneric("getTraits", function(x)standardGeneric("getTraits"))
setMethod("getTraits", "gwaswloc", function(x)
 mcols(x)$Disease.Trait)


setGeneric("subsetByChromosome", function(x, ch)standardGeneric("subsetByChromosome"))
setMethod("subsetByChromosome", "gwaswloc", function(x, ch) {
 x[ which(as.character(seqnames(x)) %in% ch) ]
})

setGeneric("subsetByTraits", function(x, tr)standardGeneric("subsetByTraits"))
setMethod("subsetByTraits", "gwaswloc", function(x, tr) {
 x[ which(getTraits(x) %in% tr) ]
})

