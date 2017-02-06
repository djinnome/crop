##jmd
##3.3.11
##collateRxnAnnots.r

collateRxnAnnots <- function(meta.ec, nc.ec, s){
    #combine nc10 data
    rxn.ec <- rbind(meta.ec, nc.ec)
    rxn.ec$nc <- as.numeric(gsub('-L2R$|-R2L$', '', rxn.ec$rxn) %in% gsub('-L2R$|-R2L$', '', nc.ec$rxn))
    #add rxns in s but not annotated
    set.diff <- setdiff(s$rxn, rxn.ec$rxn)
    add2rxn.ec <- matrix(NA, nrow=length(set.diff), ncol=ncol(rxn.ec), dimnames=list(1:length(set.diff), colnames(rxn.ec)))
    add2rxn.ec[,1] <- set.diff
    colnames(add2rxn.ec) <- colnames(rxn.ec) <- sub('React|Rxn', 'rxn', sub('Equation', 'eqn', colnames(rxn.ec)))
    add2rxn.ec[,'nc'] <- 0
    rxn.ec <- rbind(rxn.ec, add2rxn.ec)
    #delete 1 copy of rxns that are in both nc10 and meta
    rxn.ec <- rxn.ec[!duplicated(rxn.ec$rxn),]
    #delete rxns not in s
    rxn.ec <- rxn.ec[rxn.ec$rxn %in% s$rxn,]
    return(rxn.ec)
}
