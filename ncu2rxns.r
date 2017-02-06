##jmd
##1.27.12
##ncu2rxns.r
##takes vector of NCUs and outputs list of rxns

NCUvector2rxns <- function(ncu.v, gpr, gpr.name='nc10.gpr', print.nonunique.msg=FALSE){
    #set case
    ncu.v <- toupper(ncu.v)
    gpr$Genes <- toupper(gpr$Genes)
    #get rid of KOs not in GPR (accounting for rows that have multiple comma-separated NCUs)
    ncu.in.gpr <- which(apply(as.matrix(ncu.v), 1, FUN=function(x){ any(unlist(strsplit(x, split=',')) %in% gpr$Genes) }))
    cat(length(ncu.v)-length(ncu.in.gpr), 'NCUs are not in', gpr.name, 'file.', '\n')     
    ncu.v <- ncu.v[ncu.in.gpr]
    #get ko'd rxns
    ncu.rxns <- apply(as.matrix(ncu.v), 1, FUN=function(x){ ncu2rxn(ncu=unlist(strsplit(x, split=',')), gpr=gpr) })
    names(ncu.rxns) <- ncu.v
    #get rid of KOs which did not KO any rxn, but can count these as successes if our model inc. the rxn, since the success is due to GPR mapping
    ncu.no.rxn.ind <- which(sapply(ncu.rxns, FUN=function(x){ length(x)==0 }))
    cat(length(ncu.no.rxn.ind), 'NCUs are in', gpr.name, 'file but do not KO any rxns.', '\n') 
    if (length(ncu.no.rxn.ind)>=1){ ncu.rxns <- ncu.rxns[-ncu.no.rxn.ind] }
    #dups
    if (print.nonunique.msg){ cat(paste(names(ncu.rxns)[ncu.rxns %in% ncu.rxns[duplicated(ncu.rxns)]], collapse=' & '), 'have non-unique rxn sets', '\n') }
    
    return(ncu.rxns)
}

ncu2rxn <- function(ncu, gpr){
    #get rxns & enzymes associated w/ ncu
    rxns.tmp <- gpr$RXN[gpr$Genes %in% ncu]
    enzymes.tmp <- gpr$Enzyme[gpr$Genes %in% ncu]
    #get gpr w/ associated rxns which do not use KO'd enzymes
    rxns.no.ko <- gpr[(gpr$RXN %in% rxns.tmp) & !(gpr$Enzyme %in% enzymes.tmp), 'RXN']
    #KO'd rxns are setdiff
    rxns.ko <- setdiff(rxns.tmp, rxns.no.ko)
    return(rxns.ko)
}
