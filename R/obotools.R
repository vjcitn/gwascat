# obotools.R -- VJ Carey 19 March 2012
# simple models for obo term data

cleanup = function (x) 
{
#
# structure information in list read of a OBO Term recordset
# the first token in each record is regarded as a key, and the :\ is removed,
# all remaining information in the record is the value associated with the key
#
    if (!is(x, "character")) stop("expecting character vector")
    x = x[-1] # drop Term
    x = gsub(" ! .*", "", x) # deal with extraneous info in isa fields
    if (any(nchar(x) == 0)) 
        x = x[-which(nchar(x) == 0)]
    dat = strsplit(x, ": ")
    keys = sapply(dat, "[", 1)
    vals = sapply(dat, "[", -1)
    tmp = list(keys = keys, vals = vals)
    split(tmp$vals, tmp$keys)
}


obo2graphNEL = function(obo="human-phenotype-ontology.obo", 
  kill="\\[Typedef\\]") {
#
# generate a bioconductor graph graphNEL instance with
# formal term IDs as nodes, edge list defined by is_a links,
# and nodeData composed of name, def, and xref content
# 
#
 require(graph)
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
 edl = lapply(isas, function(x) list(edges=x))
 names(edl) = idvec 
 g = new("graphNEL", nodes=idvec, edgeL = edl, edgemode = "directed")
 nodeDataDefaults(g) = list(name="", def="", xref="")
 nodeData(g, nodes(g), "name") = nms
 nodeData(g, nodes(g), "def") = defs
 nodeData(g, nodes(g), "xref") = xrefs
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

