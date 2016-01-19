tt38saved.old = structure(c(957L, 699L, 649L, 323L, 294L, 249L, 248L, 245L, 222L, 
197L), .Dim = 10L, .Dimnames = structure(list(c("Obesity-related traits", 
"IgG glycosylation", "Height", "Type 2 diabetes", "Rheumatoid arthritis", 
"Crohn's disease", "Schizophrenia", "Blood metabolite levels", 
"HDL cholesterol", "Breast cancer")), .Names = ""))

# as of 18 Jan 2016

tt38saved = structure(c(957L, 822L, 699L, 340L, 294L, 249L, 248L, 245L, 222L, 
199L), .Dim = 10L, .Dimnames = structure(list(c("Obesity-related traits", 
"Height", "IgG glycosylation", "Type 2 diabetes", "Rheumatoid arthritis", 
"Crohn's disease", "Schizophrenia", "Blood metabolite levels", 
"HDL cholesterol", "Breast cancer")), .Names = ""), class = "table")

library(gwascat)
data(ebicat38)
checkTrue(identical(topTraits(ebicat38), as.table(tt38saved)))

colnames_2016jan18 = c("DATE ADDED TO CATALOG", "PUBMEDID", "FIRST AUTHOR", "DATE", 
"JOURNAL", "LINK", "STUDY", "DISEASE/TRAIT", "INITIAL SAMPLE DESCRIPTION", 
"REPLICATION SAMPLE DESCRIPTION", "REGION", "CHR_ID", "CHR_POS", 
"REPORTED GENE(S)", "MAPPED_GENE", "UPSTREAM_GENE_ID", "DOWNSTREAM_GENE_ID", 
"SNP_GENE_IDS", "UPSTREAM_GENE_DISTANCE", "DOWNSTREAM_GENE_DISTANCE", 
"STRONGEST SNP-RISK ALLELE", "SNPS", "MERGED", "SNP_ID_CURRENT", 
"CONTEXT", "INTERGENIC", "RISK ALLELE FREQUENCY", "P-VALUE", 
"PVALUE_MLOG", "P-VALUE (TEXT)", "OR or BETA", "95% CI (TEXT)", 
"PLATFORM [SNPS PASSING QC]", "CNV", "MAPPED_TRAIT", "MAPPED_TRAIT_URI"
)

