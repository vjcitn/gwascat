# obotools.R -- VJ Carey 19 March 2012
# simple models for obo term data

#iconv(df[,badinds[i]], to="ASCII", sub="*")

cleanup = function (x) 
{
#
# structure information in list read of a OBO Term recordset
# the first token in each record is regarded as a key, and the :\ is removed,
# all remaining information in the record is the value associated with the key
#
    if (!is(x, "character")) stop("expecting character vector")
    x = x[-1] # drop Term
    x = iconv(x, to="ASCII", sub="*")
    x = gsub(" ! .*", "", x) # deal with extraneous info in isa fields
    if (any(nchar(x) == 0)) 
        x = x[-which(nchar(x) == 0)]
    dat = strsplit(x, ": ")
    keys = sapply(dat, "[", 1)
    vals = sapply(dat, "[", -1)
    tmp = list(keys = keys, vals = vals)
    split(tmp$vals, tmp$keys)
}




#' convert a typical OBO text file to a graphNEL instance (using Term elements)
#' 
#' convert a typical OBO text file to a graphNEL instance (using Term elements)
#' 
#' Very rudimentary list and grep operations are used to retain sufficient
#' information to map the DAG to a graphNEL, using formal term identifiers as
#' node names and 'is-a' relationships as edges, and term names and other
#' metadata are assigned to nodeData components.
#' 
#' @aliases obo2graphNEL efo.obo.g node2uri uri2node
#' @param obo string naming a file in OBO format
#' @param kill entity types to be excluded from processing -- probably this
#' should be in a 'keep' form, but for now this works.
#' @param killTrailSp In the textual version of EFO ca. Aug 2015, there is a
#' trailing blank in the tag field defining EFO:0000001, which is not present
#' in references to this term.  Set this to TRUE to eliminate this, or graphNEL
#' construction will fail to validate.
#' @return a graphNEL instance 
#' @note The OBO for Human Disease ontology is serialized as text with this
#' package.
#' @author VJ Carey <stvjc@@channing.harvard.edu>
#' @references For use with human disease ontology,
#' \url{http://www.obofoundry.org/cgi-bin/detail.cgi?id=disease_ontology}
#' @keywords models
#' @examples
#' 
#' data(efo.obo.g)
#' requireNamespace("graph")
#' hn = graph::nodes(efo.obo.g)[1:5]
#' hn
#' graph::nodeData(efo.obo.g, hn[5])
#' 
#' @export obo2graphNEL
obo2graphNEL = function(obo="human-phenotype-ontology.obo", 
  kill="\\[Typedef\\]", killTrailSp=TRUE) {
#
# generate a bioconductor graph graphNEL instance with
# formal term IDs as nodes, edge list defined by is_a links,
# and nodeData composed of name, def, and xref content
# 
# killTrailSp introduced for EFO, EFO:0000001 has a trailing blank
# that kills the edge parsing utility
#
 if (!requireNamespace("graph")) stop("install graph package to use this function")
 lol = obo2lol(obo, kill=kill)
 xll = lapply( lol, cleanup )
 isas = lapply(xll, "[[", "is_a") #get_isas(lol)  # some can be vectors
 isas = lapply(isas, function(x) gsub(" ! .*", "", x))
#
# to retain additional content like alt_id or synonym, add code here
#
 ids = lapply(xll, "[[", "id") 
 nms = lapply(xll, "[[", "name")
 defs = lapply(xll, "[[", "def")
 xrefs = lapply(xll, "[[", "xref")
# ids = get_ids(lol)
 names(isas) = unlist(ids)
 idvec = unlist(ids)
 if (killTrailSp) idvec = gsub(" $", "", idvec)
 edl = lapply(isas, function(x) list(edges=x))
 names(edl) = idvec 
 g = new("graphNEL", nodes=idvec, edgeL = edl, edgemode = "directed")
 graph::nodeDataDefaults(g) = list(name="", def="", xref="")
 graph::nodeData(g, graph::nodes(g), "name") = nms
 graph::nodeData(g, graph::nodes(g), "def") = defs
 graph::nodeData(g, graph::nodes(g), "xref") = xrefs
 g
}

obo2lol = function(obo="human-phenotype-ontology.obo", kill="\\[Typedef\\]") {
# 
# convert the OBO file exemplified by HPO to a list of lists
# 
# note that Disease Ontology has Typedef in addition to Term elements
# and we don't want to handle these now
#
 x = readLines(obo)
 nrec = length(x)
 termstarts = grep("^\\[Term\\]$", x)
 termends = c(termstarts[-1]-1, nrec)
 tinds = lapply(1:length(termstarts), function(z) termstarts[z]:termends[z])
 termdata = lapply(tinds, function(z)x[z])
 # if kill is found in a term, then it and all subsequent records are removed
 if (nchar(kill)>0) termdata = lapply(termdata, function(x) {
    if (length(ki <- grep(kill, x))>0) {
         lx = length(x)
         x = x[-seq(ki[1],lx)]
         }
         x
    })
 termdata
}

 
# exemplar
#[Term]
#id: HP:0200053
#name: Monodactyly (feet)
#alt_id: HP:0005618
#def: "Shortening of a leg affecting only one side." [HPO:curators]
#synonym: "Asymmetric leg shortening" EXACT []
#synonym: "Asymmetric lower limb shortness" EXACT []
#xref: UMLS:C1844734 "Asymmetric lower limb shortness"
#is_a: HP:0001849 ! Oligodactyly (feet)
#created_by: sebastiankohler
#creation_date: 2011-12-02T03:41:26Z

