##jmd
##1.27.12
##merge_s.r
##merge global and organism-specific 3-column S matrices for FARM

#priority breaks ties if 2 mats have same rxn - stoichiometry could be different
#1st column must hold rxns
merge.s <- function(univ.s, org.s, col.names=c('rxn', 'compound', 'coeff'), priority=c('org', 'univ')){
    colnames(org.s) <- colnames(univ.s) <- col.names
    if (priority[1]=='org'){
        univ.dup <- univ.s[,1] %in% org.s[,1]
        s <- rbind(org.s, univ.s[!univ.dup,])
    } else {
        cat('priority set to global\n')
        org.dup <- org.s[,1] %in% univ.s[,1]
        s <- rbind(univ.s, org.s[!org.dup,])
    }
    return(s)
}

merge.smm <- function(univ.smm, org.smm, smm.id.col='FRAME'){
    #fix ncol of nc.meta st can merge
    if (length(sed <- setdiff(colnames(org.smm), colnames(univ.smm)))>0){
        for (sedi in 1:length(sed)){ univ.smm[[ sed[sedi] ]] <- NA }
        univ.smm <- univ.smm[,colnames(org.smm)]
    }
    ##combine
    smm <- rbind(org.smm, univ.smm)
    smm <- smm[!duplicated(smm[,smm.id.col]),]
    rownames(smm) <- smm[,smm.id.col]
    return(smm)
}

merge.ec <- function(univ.ec, org.ec, rxn.col='rxn', priority=c('org', 'univ')){
    if (priority[1]=='org'){ rxn.ec <- rbind(org.ec, univ.ec) } else { rxn.ec <- rbind(univ.ec, org.ec) }
    rxn.ec <- rxn.ec[!duplicated(rxn.ec[,rxn.col]),]
    rownames(rxn.ec) <- rxn.ec[,rxn.col]
    return(rxn.ec)
}
