##jmd
##9.30.11
##trace_mets.r

require('igraph')

#igraph uses zero-centering, but V(g)$name is 1-centered :(
trace.mets <- function(sp, v=NULL, from, to, pwy=NULL,
rm.mets.grep='^((A|C|G|T|U)(M|D|T)P|^NAD(P|)(H|)|FAD(H.*|)|FMN|CO-A|.+redoxin(s|)|Cytochromes.*|AMMONIA|NITR(A|I)TE|Pi|PPI|PROTON|WATER|CARBON-DIOXIDE|HYDROGEN-PEROXIDE|OXYGEN-MOLECULE)(\\[|$)' ){
    stopifnot(from %in% rownames(sp), to %in% rownames(sp), is.null(v)|all(names(v) %in% colnames(sp)))
    #subset matrix
    if (!is.null(v)){ sp <- sp[,names(v)[v>0]] }
    cat('Removing mets:', rm.mets <- sort(grep(rm.mets.grep, rownames(sp), value=TRUE)), '\n')
    sp <- sp[!(rownames(sp) %in% rm.mets),]
    nrxns <- ncol(sp); nmets <- nrow(sp)
    #bipartite mets & rxns
    nms0 <- c(rownames(sp), colnames(sp))
    sp.bip <- Matrix(0,nrow=nmets+nrxns,nmets+nrxns,dimnames=list(nms0,nms0))
    sp.bip[nmets+1:nrxns, 1:nmets] <- t(sp>0)
    sp.bip[1:nmets, nmets+1:nrxns] <- sp<0
    sp.bip[rowSums(sp.bip)!=0, rowSums(sp.bip)!=0]
    nms <- rownames(sp.bip)
    #graph
    g <- graph.adjacency(adjmatrix=sp.bip, mode='directed', diag=FALSE)
    #shortest path
    from.ind <- which(nms==from)-1; to.ind <- which(nms==to)-1
    gsp <- get.shortest.paths(graph=g, from=from.ind, to=to.ind, mode = "out")[[1]]
    gsp.nms <- nms[gsp+1]
    #structure output
    (mat <- data.frame(mets=gsp.nms[seq(1,length(gsp.nms),2)], rxns=c(gsp.nms[seq(2,length(gsp.nms),2)], '')))
    #pwys
    if (!is.null(pwy)){ mat$pwy <- pwy$Pathway.Name[match(mat$rxns, pwy$rxn)] }
    return(mat)
}
