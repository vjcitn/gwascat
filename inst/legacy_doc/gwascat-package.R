

#' 
#' representing and modeling data in the NHGRI GWAS catalog, using GRanges and
#' allied infrastructure
#' 
#' \tabular{ll}{ Package: \tab gwascat\cr Version: \tab 1.7.3\cr Suggests: \tab
#' \cr Depends: \tab R (>= 3.0.0), methods, IRanges, GenomicRanges\cr Imports:
#' \tab \cr License: \tab Artistic-2.0\cr LazyLoad: \tab yes\cr }
#' 
#' Index: \preformatted{ gwaswloc-class Class '"gwaswloc"' }
#' 
#' The GWAS catalog management has migrated to EMBL/EBI.  Use data(ebicat38)
#' for an image dated 3 August 2015.  Use makeCurrentGwascat() to get a more
#' recent image.  Use data(ebicat37) for a GRCh37 (or hg19) liftOver result.
#' Use data(ebicat37UCSC) for an image with \code{hg19} as genome tag and UCSC
#' seqnames.
#' 
#' The data objects
#' 
#' 'g17SM' 'gg17N' 'gw6.rs_17' 'low17' 'rules_6.0_1kg_17' 'gwrngs'
#' 
#' are described in vignettes.
#' 
#' The DataFrame function is imported from IRanges.
#' 
#' The \code{\link[snpStats]{SnpMatrix-class}} is used to represent data
#' related to rule-based imputation, using the
#' \code{\link[snpStats]{impute.snps}} function.
#' 
#' si.hs.38 is a \code{\link[Seqinfo]{Seqinfo-class}} instance for hg38.
#' 
#' @name gwascat-package
#' @aliases gwascat-package gwascat g17SM gg17N gw6.rs_17 low17
#' rules_6.0_1kg_17 DataFrame gwrngs19 gwrngs38 ebicat38 impute.snps
#' SnpMatrix-class si.hs.38
#' @docType package
#' @author VJ Carey <stvjc@@channing.harvard.edu>
#' 
#' Maintainer: VJ Carey <stvjc@@channing.harvard.edu>
#' @references \url{http://www.ebi.ac.uk/gwas/}
#' 
#' Partial support from the Computational Biology Group at Genentech, Inc.
#' @keywords package
#' @examples
#' 
#'  data(ebicat38)
#'  ebicat38
#' 
NULL





#' 
#' data on 1000 genomes SNPs that 'tag' GWAS catalog entries
#' 
#' This GRanges instance includes locations for 297000 1000 genomes SNP that
#' have been identified as exhibiting LD with NHGRI GWAS SNP as of September
#' 2013.  The tagid field tells the name of the tagging SNP, the baseid field
#' is the SNP identifier for the GWAS catalog entry, the R2 field tells the
#' value of R-squared relating the distributions of the tagging SNP and the
#' GWAS entry.  Only tagging SNP with R-squared 0.8 or greater are included. A
#' self-contained R-based procedure should emerge in 2014.
#' 
#' @name gwastagger
#' @docType data
#' @format The format is:\cr Formal class 'GRanges' [package "GenomicRanges"]
#' with 6 slots\cr ..@ seqnames :Formal class 'Rle' [package "IRanges"] with 4
#' slots\cr .. .. ..@ values : Factor w/ 24 levels "chr1","chr2",..: 1 2 3 4 5
#' 6 7 8 9 10 ...\cr .. .. ..@ lengths : int [1:22] 24042 23740 21522 14258
#' 14972 34101 12330 11400 8680 15429 ...\cr .. .. ..@ elementMetadata: NULL\cr
#' .. .. ..@ metadata : list()\cr ..@ ranges :Formal class 'IRanges' [package
#' "IRanges"] with 6 slots\cr .. .. ..@ start : int [1:297579] 986111 988364
#' 992250 992402 995669 999686 1005579 1007450 1011209 1011446 ...\cr .. .. ..@
#' width : int [1:297579] 1 1 1 1 1 1 1 1 1 1 ...\cr .. .. ..@ NAMES : NULL\cr
#' .. .. ..@ elementType : chr "integer"\cr .. .. ..@ elementMetadata: NULL\cr
#' .. .. ..@ metadata : list()\cr ..@ strand :Formal class 'Rle' [package
#' "IRanges"] with 4 slots\cr .. .. ..@ values : Factor w/ 3 levels
#' "+","-","*": 3\cr .. .. ..@ lengths : int 297579\cr .. .. ..@
#' elementMetadata: NULL\cr .. .. ..@ metadata : list()\cr ..@
#' elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots\cr
#' .. .. ..@ rownames : NULL\cr .. .. ..@ nrows : int 297579\cr .. .. ..@
#' listData :List of 3\cr .. .. .. ..$ tagid : chr [1:297579] "rs28479311"
#' "rs3813193" "chr1:992250" "rs60442576" ...\cr .. .. .. ..$ R2 : num
#' [1:297579] 0.938 0.994 0.969 1 1 ...\cr .. .. .. ..$ baseid: chr [1:297579]
#' "rs3934834" "rs3934834" "rs3934834" "rs3934834" ...\cr .. .. ..@ elementType
#' : chr "ANY"\cr .. .. ..@ elementMetadata: NULL\cr .. .. ..@ metadata :
#' list()\cr ..@ seqinfo :Formal class 'Seqinfo' [package "GenomicRanges"] with
#' 4 slots\cr .. .. ..@ seqnames : chr [1:24] "chr1" "chr2" "chr3" "chr4"
#' ...\cr .. .. ..@ seqlengths : int [1:24] 249250621 243199373 198022430
#' 191154276 180915260 171115067 159138663 146364022 141213431 135534747 ...\cr
#' .. .. ..@ is_circular: logi [1:24] FALSE FALSE FALSE FALSE FALSE FALSE
#' ...\cr .. .. ..@ genome : chr [1:24] "hg19" "hg19" "hg19" "hg19" ...\cr ..@
#' metadata : list()\cr
#' @source NHGRI GWAS catalog; plink is used with the 1000 genomes VCF in a
#' perl routine by Michael McGeachie, Harvard Medical School;
#' @keywords datasets
#' @examples
#' 
#' data(gwastagger)
#' gwastagger[1:5]
#' data(ebicat37)
#' mean(ebicat37$SNPS %in% gwastagger$baseid)
#' # ideally, all GWAS SNP would be in our tagging ranges as baseid
#' query <- setdiff(ebicat37$SNPS, gwastagger$baseid)
#' # relatively recent catalog additions
#' sort(table(ebicat37[which(ebicat37$SNPS %in% query)]$DATE.ADDED.TO.CATALOG), decreasing=TRUE)[1:10]
#' 
NULL





#' Class \code{"gwaswloc"}
#' 
#' A container for GRanges instances representing information in the NHGRI GWAS
#' catalog.
#' 
#' @name gwaswloc-class
#' @importFrom readr read_tsv
#' @importFrom S4Vectors DataFrame elementNROWS queryHits subjectHits values "values<-" metadata "metadata<-"
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges findOverlaps overlapsAny 
#' @importFrom Rsamtools TabixFile
#' @importFrom ggplot2 aes
#' @importFrom rtracklayer import
#' @importFrom methods as callNextMethod is
#' @importFrom utils data read.delim
#' @importFrom Seqinfo seqlevels seqlengths "seqlengths<-" "seqlevels<-" "seqinfo" "seqinfo<-" seqnames
#' @importFrom GenomeInfoDb "seqlevelsStyle<-"
#' @aliases gwaswloc-class [,gwaswloc,ANY,ANY,ANY-method [,gwaswloc,ANY-method
#' [,gwaswloc-method show,gwaswloc-method getRsids,gwaswloc-method
#' getTraits,gwaswloc-method subsetByChromosome,gwaswloc-method
#' subsetByTraits,gwaswloc-method getTraits getRsids subsetByChromosome
#' subsetByTraits
#' @docType class
#' @note In gwascat package 1.9.6 and earlier, the globally accessible
#' \code{gwaswloc} instance \code{gwrngs} was created upon attachment.  This is
#' no longer the case.
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("gwaswloc", ...)}. Any GRanges instance can be supplied.
#' @author VJ Carey <stvjc@@channing.harvard.edu>
#' @references %% ~~put references to the literature/web site here~~
#' \url{http://www.ebi.ac.uk/gwas/}
#' @keywords classes
#' @examples
#' 
#' showClass("gwaswloc")
#' 
NULL





#' internal data frame for NHGRI GWAS catalog
#' 
#' convenience container for imported table from NHGRI GWAS catalog
#' 
#' 
#' @name gwdf_2012_02_02
#' @aliases gwdf_2014_09_08
#' @docType data
#' @format A data frame with 17832 observations on the following 34 variables.
#' @return a DataFrame with (character) columns
#' "Date Added to Catalog",
#' "PUBMEDID",
#' "First Author",
#' "Date",
#' "Journal",
#' "Link",
#' "Study",
#' "Disease/Trait",
#' "Initial Sample Size",
#' "Replication Sample Size",
#' "Region",
#' "Chr_id",
#' "Chr_pos",
#' "Reported Gene(s)",
#' "Mapped_gene",
#' "Upstream_gene_id",
#' "Downstream_gene_id",
#' "Snp_gene_ids",
#' "Upstream_gene_distance",
#' "Downstream_gene_distance",
#' "Strongest SNP-Risk Allele",
#' "SNPs",
#' "Merged",
#' "Snp_id_current",
#' "Context",
#' "Intergenic",
#' "Risk Allele Frequency",
#' "p-Value",
#' "Pvalue_mlog",
#' "p-Value (text)",
#' "OR or beta",
#' "95% CI (text)",
#' "Platform \[SNPs passing QC\]",
#' "CNV"
#' @note In versions prior to 1.9.6, The .onAttach function specifies which
#' data frame is transformed to GRanges.  This is now managed manually.
#' @source \url{http://www.ebi.ac.uk/gwas/}
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' data(gwdf_2014_09_08)
#' # try gwascat:::gwdf2GRanges on this data.frame
#' }
#' 
NULL





#' location information for 10000 SNPs probed on Affy GW 6.0
#' 
#' location information for 10000 SNPs probed on Affy GW 6.0
#' 
#' extracted from pd.genomewidesnp.6 v 1.4.0; for demonstration purposes
#' 
#' @name locon6
#' @docType data
#' @format A data frame with 10000 observations on the following 3 variables.
#' \describe{ \item{dbsnp_rs_id}{a character vector}
#' \item{chrom}{a character vector} \item{physical_pos}{a
#' numeric vector} }
#' @keywords datasets
#' @examples
#' 
#' data(locon6)
#' str(locon6) 
#' 
NULL
