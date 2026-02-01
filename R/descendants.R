
descendants = function(g, root, depth=1) {
 #
 # return vector of nodes at distance at most depth from root
 # counting edges (not weights)
 # include root
 #
 if (!requireNamespace("graph")) stop("install graph package to use this function")
 alln = list()
 alln[[1]] = root
 for (i in 1:depth) 
   alln[[i+1]] = unlist(graph::adj(g, unlist(alln)))
 unique(unlist(alln))
}

parentingDF = function(gnel, nodes, useNodeData=NULL) {
 if (!requireNamespace("graph")) stop("install graph package to use this function")
 if (!is(gnel, "graphNEL")) stop("currently for graphNEL input")
 ng = graph::subGraph(nodes, gnel)
 mapperf = function(x)x
 mapper = NULL
 if (!is.null(useNodeData)) {
   if (length(useNodeData)>1 | !is.character(useNodeData)) stop("useNodeData must be atomic character")
   nd = graph::nodeData(ng, , useNodeData)
   mapper = unlist(nd)
   names(mapper) = names(nd)
   mapperf = function(x) mapper[x]
 }
 ng = as(ng, "graphBAM")
 ftmat = graph::extractFromTo(ng)  # data.frame
 ftmat = lapply(ftmat[1:2], as.character)
 ftmat = do.call(cbind, ftmat)
 formal = apply(ftmat, 1, function(x) paste(x, collapse="_"))
 colnames(ftmat)[1] = "Parent"
 dd = dim(ftmat)
 dn = dimnames(ftmat)
 x = mapperf( as.character(ftmat) )
 mx = matrix(x, nrow=dd[1])
 dimnames(mx) = dn
 data.frame(cbind(mx, formal))
}
 
