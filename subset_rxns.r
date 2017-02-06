##jmd
##8.16.11
##subset_rxns.r

#don't remove 'ARYLFORMAMIDASE-RXN-L2R' anymore

subset.rxns <- function(rxns, filt='Neurospora/nc10.filtered', wrong.dir='nc_wrong_dir_rxns.txt', meta.add='meta_add_rxns.txt', filt.add='filt_add_rxns.txt',
wrong.dir.add='wrong_dir_add_rxns.txt',
rid.rxns=c('SUCROSE-SYNTHASE-RXN-CPD-12575/BETA-D-FRUCTOSE//SUCROSE/UDP/PROTON.46.-R2L', 'XANPRIBOSYLTRAN-RXN-R2L',
'FRUCTOSE-6-PHOSPHATE-PHOSPHOKETOLASE-RXN-L2R', 'GLUCONOKIN-RXN-L2R', 'CIT-TRANS-RXN-L2R')){

    ##read
    filt.rxns <- read.delim(filt)
    wrong.dir.rxns <- read.delim(wrong.dir)
    #add rxns
    meta.add.rxns <- read.delim(meta.add)
    filt.add.rxns <- read.delim(filt.add)
    wrong.dir.add.rxns <- read.delim(wrong.dir.add)
    
    ##set rxns
    ss0 <- setdiff(rownames(rxns)[rxns$nc==1], c(filt.rxns$rxn, wrong.dir.rxns$rxn))
    ss1 <- c(ss0, rownames(meta.add.rxns), rownames(filt.add.rxns), rownames(wrong.dir.add.rxns))
    ss2 <- setdiff(ss1, rid.rxns)
    ss <- unique(ss2[!is.na(ss2)])
    
    return(ss)
}
