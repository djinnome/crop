##jmd
##10.18.11
##get_biomass_mass.r

get.bm.mass <- function(bm, smm){
    #don't include mass of atp & fadh2 in left-side of biomass compositions
    bm[bm$comp %in% c('ATP', 'FADH2') & bm$rxn %in% c('biomass', 'FullBiomassComposition'), 'coeff'] <- 0
    wts <- unlist(apply(as.matrix(unique(bm$rxn)), 1, FUN=function(set){ 
        bm2 <- bm[bm$rxn %in% set & bm$coeff<0,]
        mets <- gsub('^Charged-|-tRNAs$', '', bm2$comp) 
        mass <- smm[mets, 'MASS']/1000
        #if mass is NA, assume it's 1, since we don't have any instances in biomass composition that aren't in smm
        mass[is.na(mass)] <- 1
        ret <- -bm2$coeff %*% mass
    })) #end apply fcn
    names(wts) <- unique(bm$rxn)
    return(wts)
}
