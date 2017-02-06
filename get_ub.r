##jmd
##9.13.11
##get_ub.r
##code to get tight bounds using fcns
#called from get_ub_script.r

get.ub <- function(a, eps=10**-5, ex.rxns, sense='E'){
    #instantiate
    lb <- rep(0, ncol(a)); ub <- rep(Inf, ncol(a)); obj <- rep(-eps, ncol(a))
    names(obj) <- names(lb) <- names(ub) <- colnames(a)
    #lb['biomass'] <- 1
    ub[names(dub.nona)] <- 1/dub.nona + 0.01
    #want this flux to be used for biomass, so don't allow flux thru export rxns, but this doesn't prevent mass-depleting rxns
    #ub[ex.rxns] <- 0
    
    ub.tmp <- ub
    #ret <- rep(NA, ncol(a)); names(ret) <- colnames(a)
    #for (i in 201:300){
        obj[i] <- 1
        #want ub[i] lower than other ub's, since don't want route from rxn to biomass to inc. mult paths due to insufficient flux along path
        ub.tmp[i] <- 10**5
        max.v <- Rcplex(cvec=obj, Amat=a, bvec=rep(0, nrow(a)), lb=lb, ub=ub.tmp, sense=sense, objsense='max', control=list(trace=0))
        names(max.v$xopt) <- colnames(a); v <- max.v$xopt
        cat('CPLEX status:', max.v$status, 'w/ max objective:', v[i], '\n')
        
        c2 <- obj; c2[c2<0] <- 0
        u2 <- v; u2[ex.rxns] <- 10**6
        min.v <- Rcplex(cvec=c2, Amat=a, bvec=rep(0, nrow(a)), lb=lb, ub=u2, sense=sense, objsense='min', control=list(trace=0))
        names(min.v$xopt) <- colnames(a); v2 <- min.v$xopt
        cat('CPLEX status:', min.v$status, 'w/ min objective:', min.v$obj, '\n \n')
        if (min.v$status==1) ret[i] <- min.v$obj
        
        ub.tmp <- ub; obj <- rep(-eps, ncol(a))
    }#end for
    return(ret)
}
