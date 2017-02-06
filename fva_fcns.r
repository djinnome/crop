##jmd
##8.5.11
##fva_fcns.r
#use goal programming for faster maximization, to see if rxn can carry flux or not
#assumes sv>=0 & rev rxns split into 2 irrev rxns

require('Matrix')
require('Rcplex')

max.flux <- function(a, sense='E', rxns.ss=colnames(a), v.ub=1, ub=rep(10^4, ncol(a)), allow.all.trans=TRUE){
    stopifnot(rxns.ss %in% colnames(a))
    flux <- lb <- obj <- numeric(ncol(a))
    names(lb) <- names(ub) <- names(flux) <- names(obj) <- colnames(a)
    ub.tmp <- ub; obj.tmp <- obj
    if (allow.all.trans){ a <- a[-grep('\\[CCO-EXTRACELLULAR\\]$', rownames(a)),] }
    for (rxn.tmp in rxns.ss){
        obj.tmp[rxn.tmp] <- 1
        ub.tmp[rxn.tmp] <- v.ub
        max.rxn <- Rcplex(cvec=obj.tmp, Amat=a, lb=lb, ub=ub.tmp, bvec=numeric(nrow(a)), objsense='max', sense=sense, control=list(trace=0))
        if (max.rxn$stat==1) flux[rxn.tmp] <- max.rxn$obj else flux[rxn.tmp] <- NA
        #reset
        ub.tmp <- ub; obj.tmp <- obj
    }
    return(flux[rxns.ss])
}

#use ub=Inf, the Rcplex default
min.flux <- function(a, fba.obj.rxn='biomass', min.bm=50, sense='E', rxns.min=colnames(a)){
    stopifnot(c(fba.obj.rxn, rxns.min) %in% colnames(a))
    lb <- obj <- flux <- numeric(ncol(a))
    names(lb) <- names(obj) <- names(flux) <- colnames(a)
    lb[fba.obj.rxn] <- min.bm
    #trans.rxns <- grep('^TRANS-RXN', colnames(a), value=TRUE)
    #ub[trans.rxns] <- 10**7
    obj.tmp <- obj
    for (i in rxns.min){
        obj.tmp[i] <- 1
        min.rxn <- Rcplex(cvec=obj.tmp, Amat=a, lb=lb, bvec=numeric(nrow(a)), sense=sense, objsense='min', control=list(trace=0))
        if (min.rxn$stat==1) flux[i] <- min.rxn$obj else flux[i] <- NA
        #reset
        obj.tmp <- obj
    }
    return(flux[rxns.min])
}

##on 9.6.11, this gave same number of rxns w/ no flux as max.flux in fva_fcns.r
#CPU system.time comparison: max.n.flux=0.004; genome-scale max.flux=6.5. % savings>99.9%
#max slack var t st S_intra*v=0, v>=t, 1>=t>=0.
max.n.flux <- function(sp, v.ub=Inf, se='E', slack.goal.ub=1, allow.all.trans=TRUE){
    fba.se <- character(nrow(sp)); names(fba.se) <- rownames(sp)
    fba.se[1:nrow(sp)] <- se
    #if allow reversible transport of all extracellular compounds
    if (allow.all.trans){ 
        # this will exclude purely extracell rxns, eg invertase
        sp <- sp[-grep('\\[CCO-EXTRACELLULAR\\]$', rownames(sp)),] 
        fba.se <- fba.se[rownames(sp)]
    }
    nrxns <- ncol(sp); nmets <- nrow(sp)
    #decision vars are: v slack
    Amat <- rBind(cBind(sp, Matrix(0, nrow=nmets, ncol=nrxns)),
            cBind(Diagonal(n=nrxns), -1*Diagonal(n=nrxns)))
    ub <- rep(c(v.ub, slack.goal.ub), each=nrxns)
    obj <- rep(c(0,1), each=nrxns)
    sense <- c(fba.se, rep('G', nrxns))
    
    max.nv <- Rcplex(Amat=Amat, cvec=obj, bvec=numeric(nmets+nrxns), sense=sense, ub=ub, lb=numeric(2*nrxns), objsense='max')
    names(max.nv$xopt) <- c(colnames(sp), paste('slack', colnames(sp), sep='_'))
    cat('CPLEX status:', max.nv$status, 'w/ n no flux rxns:', sum(max.nv$xopt[-(1:nrxns)]<10**-6), '\n')
    #some fluxes are epsilon close to goal
    v <- max.nv$xopt[1:nrxns]; v[v>=slack.goal.ub-10**-6] <- 1
    return(v)
}
