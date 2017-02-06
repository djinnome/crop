##jmd
##8.11.11
##minFlux4met.r
##min fluxes to accumulate 1 unit of a (set of) metabolites

minFlux4met <- function(a, met, ko.rxns=NULL, sense='E', ub=rep(10^3, ncol(a)), ...){
    stopifnot(met %in% rownames(a))
    
    lb <- rep(0, ncol(a)); obj <- rep(1, ncol(a))
    names(obj) <- names(lb) <- names(ub) <- colnames(a)
    ub[ko.rxns] <- 0
    
    rhs <- numeric(nrow(a)); names(rhs) <- rownames(a); rhs[met] <- 1
    
    min.v <- Rcplex(cvec=obj, Amat=a, bvec=rhs, ub=ub, lb=lb, sense=sense, objsense='min', ...)
    names(min.v$xopt) <- colnames(a)
    
    cat('CPLEX status:', min.v$status, '\n')
    return(min.v)
}#end fcn
