##jmd
##2.27.12
##run_get_md_constraints_3c.r


run.get.md.constraints.3c <- function(sp, mets.intra, md.f.ub){
    #make sp.intra, which has same number of rxns as sp
    sp.intra <- sp[mets.intra,]
    sp.intra.norm <- abs(sp.intra)/rowSums(abs(sp.intra))        
    sp.intra.3c <- dense2sp(sp.intra)
    sp.intra.norm.3c <- dense2sp(sp.intra.norm)
    
    stopifnot(ncol(sp)==ncol(sp.intra), dim(sp.intra)==dim(sp.intra.norm)) 
    
    ret <- get.md.constraints.3c(sp.intra.3c=sp.intra.3c, sp.intra.norm.3c=sp.intra.norm.3c, rnames.sp.intra=rownames(sp.intra), cnames.sp.intra=colnames(sp.intra),
    md.f.ub=md.f.ub)
    return(ret)
}
