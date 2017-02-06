##jmd
##9.25.11
##check_new_media.r

check.new.media <- function(sp, ub=named.vec(10**3, names=colnames(sp)), supp.ub=20, rm.nuts=NULL, add.nuts=NULL, ko=NULL, ctrl=list(trace=0), ...){
    sp2 <- changeNuts(sp, rm.nuts=rm.nuts, add.nuts=add.nuts)
    
    if (!is.null(ko)){
        #if(!all(ko %in% colnames(sp))){ cat('Not all KOs in matrix \n') }
        ko <- intersect(ko, colnames(sp))
    }
    
    ub2 <- named.vec(10**3, names=colnames(sp2))
    ub2[intersect(names(ub2), names(ub))] <- ub[intersect(names(ub2), names(ub))]
    ub2[setdiff(names(ub2), names(ub))] <- supp.ub
        
    f <- fba.na(a=sp2, fba.ub=ub2, ko=ko, control=ctrl, ...)
    return(f)
}
