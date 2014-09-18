\name{gwastagger}
\alias{gwastagger}
\docType{data}
\title{
data on 1000 genomes SNPs that 'tag' GWAS catalog entries
}
\description{
data on 1000 genomes SNPs that 'tag' GWAS catalog entries
}
\usage{data(gwastagger)}
\format{
  The format is:\cr
Formal class 'GRanges' [package "GenomicRanges"] with 6 slots\cr
  ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 4 slots\cr
  .. .. ..@ values         : Factor w/ 24 levels "chr1","chr2",..: 1 2 3 4 5 6 7 8 9 10 ...\cr
  .. .. ..@ lengths        : int [1:22] 24042 23740 21522 14258 14972 34101 12330 11400 8680 15429 ...\cr
  .. .. ..@ elementMetadata: NULL\cr
  .. .. ..@ metadata       : list()\cr
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots\cr
  .. .. ..@ start          : int [1:297579] 986111 988364 992250 992402 995669 999686 1005579 1007450 1011209 1011446 ...\cr
  .. .. ..@ width          : int [1:297579] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. .. ..@ NAMES          : NULL\cr
  .. .. ..@ elementType    : chr "integer"\cr
  .. .. ..@ elementMetadata: NULL\cr
  .. .. ..@ metadata       : list()\cr
  ..@ strand         :Formal class 'Rle' [package "IRanges"] with 4 slots\cr
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 3\cr
  .. .. ..@ lengths        : int 297579\cr
  .. .. ..@ elementMetadata: NULL\cr
  .. .. ..@ metadata       : list()\cr
  ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots\cr
  .. .. ..@ rownames       : NULL\cr
  .. .. ..@ nrows          : int 297579\cr
  .. .. ..@ listData       :List of 3\cr
  .. .. .. ..$ tagid : chr [1:297579] "rs28479311" "rs3813193" "chr1:992250" "rs60442576" ...\cr
  .. .. .. ..$ R2    : num [1:297579] 0.938 0.994 0.969 1 1 ...\cr
  .. .. .. ..$ baseid: chr [1:297579] "rs3934834" "rs3934834" "rs3934834" "rs3934834" ...\cr
  .. .. ..@ elementType    : chr "ANY"\cr
  .. .. ..@ elementMetadata: NULL\cr
  .. .. ..@ metadata       : list()\cr
  ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomicRanges"] with 4 slots\cr
  .. .. ..@ seqnames   : chr [1:24] "chr1" "chr2" "chr3" "chr4" ...\cr
  .. .. ..@ seqlengths : int [1:24] 249250621 243199373 198022430 191154276 180915260 171115067 159138663 146364022 141213431 135534747 ...\cr
  .. .. ..@ is_circular: logi [1:24] FALSE FALSE FALSE FALSE FALSE FALSE ...\cr
  .. .. ..@ genome     : chr [1:24] "hg19" "hg19" "hg19" "hg19" ...\cr
  ..@ metadata       : list()\cr
}
\details{
This GRanges instance includes locations for 297000 1000 genomes SNP that have been identified
as exhibiting LD with NHGRI GWAS SNP as of September 2013.  The tagid field tells the name
of the tagging SNP, the baseid field is the SNP identifier for the GWAS catalog entry,
the R2 field tells the value of R-squared relating the distributions of the tagging SNP
and the GWAS entry.  Only tagging SNP with R-squared 0.8 or greater are included.
A self-contained R-based procedure should emerge in 2014.
}
\source{
NHGRI GWAS catalog; plink is used with the 1000 genomes VCF in a perl routine
by Michael McGeachie, Harvard Medical School; 
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(gwastagger)
gwastagger[1:5]
data(gwrngs19)
mean(gwrngs19$SNPs \%in\% gwastagger$baseid)
# ideally, all GWAS SNP would be in our tagging ranges as baseid
query <- setdiff(gwrngs19$SNPs, gwastagger$baseid)
# relatively recent catalog additions
sort(table(gwrngs19[which(gwrngs19$SNPs \%in\% query)]$Date.Added.to.Catalog), decreasing=TRUE)[1:10]
}
\keyword{datasets}
