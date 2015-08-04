tt38saved = structure(c(957L, 699L, 649L, 323L, 294L, 249L, 248L, 245L, 222L, 
197L), .Dim = 10L, .Dimnames = structure(list(c("Obesity-related traits", 
"IgG glycosylation", "Height", "Type 2 diabetes", "Rheumatoid arthritis", 
"Crohn's disease", "Schizophrenia", "Blood metabolite levels", 
"HDL cholesterol", "Breast cancer")), .Names = ""))

library(gwascat)
data(ebicat38)
checkTrue(identical(topTraits(ebicat38), tt38saved))
