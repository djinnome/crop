##jmd
##9.22.11
##check_rxn.r

check.rxn.inputs <- function(sp, rxn.tmp){
    ssm <- ss.mat(sp, c=rxn.tmp)
    f <- FBA(sp, goals=names(ssm)[ssm<0])
    return(f$x[paste(names(ssm)[ssm<0], 'GOAL-RXN-L2R', sep='-')])
}

check.rxn.outputs <- function(sp, rxn.tmp){
    ssm <- ss.mat(sp, c=rxn.tmp)
    outs <- names(ssm)[ssm>0]
    outs.trans <- paste(outs, 'TRANS-RXN-L2R', sep='-')
    
    #to test if output mets can be consumed, 1. add importers for them & 2. check if these can carry flux
    new.rxns <- matrix(0, nrow=nrow(sp), ncol=length(outs), dimnames=list(rownames(sp), outs.trans))
    new.rxns[cbind(outs, outs.trans)] <- 1
    
    sp2 <- cBind(sp, Matrix(new.rxns))
    f <- FBA(sp2, goals=outs.trans)
    return(f$x[outs.trans])
}
