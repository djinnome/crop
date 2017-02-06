##jmd
##7.18.2011
##fluxes_of_met.r

#flux vector v must have names representing rxns that match s.sp
flux.of.met <- function(a, v, met){
    stopifnot(met %in% rownames(a), names(v)==colnames(a))
    v <- abs(v[v!=0])
    met.rxns <- colnames(a)[a[met,]!=0]
    rxns.oi <- intersect(names(v), met.rxns)
    v.out <- a[met, rxns.oi] * v[rxns.oi]
    return(v.out)
}
