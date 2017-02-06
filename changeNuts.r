##jmd
##11.16.11
##changeNuts.r

changeNuts <- function(sp, rm.nuts=NULL, add.nuts=NULL){
    if (!is.null(rm.nuts)){ 
        rm.nuts <- gsub('\\[CCO-.+\\]$', '', rm.nuts)
        rm.trans.rxns <- paste(rm.nuts, 'TRANS-RXN-L2R', sep='-')
        if (!all(rm.trans.rxns %in% colnames(sp))){
            cat('rm.nuts not in sp are:', setdiff(rm.trans.rxns, colnames(sp)), '\n') 
            rm.trans.rxns <- intersect(rm.trans.rxns, colnames(sp))
        }
        sp.cn <- sp[,!(colnames(sp) %in% rm.trans.rxns)]
    } else {
        sp.cn <- sp
    }#end else
    
    if (!is.null(add.nuts)){ 
        #add.nuts <- gsub('\\[CCO-.+\\]$', '', add.nuts)
        if(!all(add.nuts %in% rownames(sp))){ 
            cat('add.nuts not in sp are:', setdiff(add.nuts, rownames(sp)), '\n') 
            add.nuts <- intersect(add.nuts, rownames(sp))
        }
        add.trans.rxns <- paste(gsub('\\[CCO-.+\\]$', '', add.nuts), 'TRANS-RXN-L2R', sep='-')
        new.mat <- matrix(0, nrow=nrow(sp), ncol=length(add.trans.rxns), dimnames=list(rownames(sp), add.trans.rxns))
        new.mat[cbind(add.nuts, add.trans.rxns)] <- 1
        sp.cn <- cBind(sp.cn, Matrix(new.mat))
    }
    return(sp.cn)
}
