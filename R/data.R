#  ‘ebicat_2020_04_30’ ‘g17SM’ ‘gg17N’ ‘gw6.rs_17’ ‘gwastagger’ ‘locon6’
#  ‘low17’ ‘si.hs.37’ ‘si.hs.38’

#' serialized gwaswloc instance from april 30 2020, sample of 50000 records
#' @usage data(ebicat_2020_04_30)
#' @docType data
#' @format gwaswloc instance
"ebicat_2020_04_30"

#' SnpMatrix instance from chr17
#' @docType data
#' @usage data(g17SM)
#' @format snpStats SnpMatrix instance
"g17SM"

#' genotype matrix from chr17 1000 genomes
#' @docType data
#' @usage data(gg17N)
#' @format matrix
#' @examples
#' data(gg17N)
#' gg17N[1:4,1:4]
"gg17N"

#' character vector of rs numbers for SNP on chr17
#' @docType data
#' @usage data(gw6.rs_17)
#' @format character vector
"gw6.rs_17"

#' GRanges with LD information on 9998 SNP
#' @docType data
#' @usage data(gwastagger)
#' @format GRanges
"gwastagger"

#' location data for 10000 SNP
#' @usage data(locon6)
#' @docType data
#' @format data.frame, coordinates are hg19
"locon6"

#' SnpMatrix instance from chr17
#' @docType data
#' @usage data(low17)
#' @format snpStats SnpMatrix instance
"low17"

#' Seqinfo for GRCh37
#' @docType data
#' @usage data(si.hs.37)
#' @format Seqinfo instance
"si.hs.37" 

#' Seqinfo for GRCh38
#' @docType data
#' @usage data(si.hs.38)
#' @format Seqinfo instance
"si.hs.38" 

#' image of locon6 in GRanges, lifted over to hg38
#' @docType data
#' @usage data(gr6.0_hg38)
#' @format GRanges instance
"gr6.0_hg38"
