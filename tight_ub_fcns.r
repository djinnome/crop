##jmd
##9.2.11
##tight_ub_fcns.r

#dual.ub.milp treats beta of rxn.tmp as 0, so n.ko includes rxn.tmp, whereas in dum.cut beta of rxn.tmp = 1
#because dum.cut uses cut.met.mat, setting n.ko in dum.cut may cause infeasibility

require('Matrix')
require('Rcplex')

#lp for only nec.rxns
dual.ub.nec <- function(sp, nec.rxns){
    dub <- rep(NA, length(nec.rxns)); names(dub) <- nec.rxns
    for (i in 1:length(nec.rxns)){
        ub <- obj <- numeric(ncol(sp)+nrow(sp)); lb <- rep(-Inf, ncol(sp)+nrow(sp))
        names(obj) <- names(ub) <- names(lb) <- paste(rep(c('mu','l'), times=c(nrow(sp), ncol(sp))), c(rownames(sp), colnames(sp)), sep='_')
        obj[paste('l', nec.rxns[i], sep='_')] <- 1
        ub[paste('l', nec.rxns[i], sep='_')] <- Inf
        #if (sense=='G') lb[-(1:nrow(sp))] <- -Inf
        #mu <= Inf
        ub[1:nrow(sp)] <- Inf
        #l.bm=0
        lb['l_biomass'] <- ub['l_biomass'] <- 0
        
        a <- cBind(t(sp), -1*Diagonal(ncol(sp)))
        dimnames(a) <- list(colnames(sp), names(obj))
        fba.obj <- numeric(ncol(sp)); names(fba.obj) <- colnames(sp); fba.obj['biomass'] <- 1
        max.dub <- Rcplex(cvec=obj, Amat=a, bvec=-fba.obj, sense='E', objsense='min', ub=ub, lb=lb, control=list(trace=0))
        if (max.dub$stat==1) dub[i] <- max.dub$obj
    }
    return(dub)
}

#bilevel milp
dual.ub.milp <- function(a, rxns.ss=colnames(a), n.ko=1, ctrl=list(tilim=5), fba.ub=rep(10**5, ncol(a)), lam.ub=rep(10**3, ncol(a))){

    nrxns <- ncol(a); nmets <- nrow(a)
    fba.obj <- Matrix(0, nrow=1, ncol=nrxns, dimnames=list(1, colnames(a))); fba.obj[1,'biomass'] <- 1
    names(fba.ub) <- names(lam.ub) <- colnames(a); trans.rxns <- grep('^TRANS-RXN', colnames(a), value=TRUE); fba.ub[trans.rxns] <- 10**3 * mean(fba.ub)
    
    #decision vars are: v mu lambda beta
    Amat <- rBind(cBind(Diagonal(x=numeric(nrxns)), t(a), -1*Diagonal(n=nrxns), Diagonal(x=numeric(nrxns))),
    cBind(Matrix(0, nrow=nrxns, ncol=nrxns+nmets), Diagonal(n=nrxns), Diagonal(x=lam.ub)),
    cBind(a, Matrix(0, nrow=nmets, ncol=2*nrxns+nmets)),
    cBind(Diagonal(n=nrxns), Matrix(0, nrow=nrxns, ncol=nmets+nrxns), Diagonal(x=-fba.ub)),
    cBind(-fba.obj, Matrix(0, nrow=1, ncol=2*nrxns+nmets)),
    cBind(Matrix(rep(c(0, 1), times=c(2*nrxns+nmets, nrxns)), nrow=1, ncol=3*nrxns+nmets)))
    
    #names
    colnames(Amat) <- c(colnames(a), paste('mu', rownames(a), sep='_'), paste('l', colnames(a), sep='_'), paste('beta', colnames(a), sep='_'))
    rownames(Amat) <- c(paste('dual_eq', colnames(a), sep='_'), paste('l_bounds', colnames(a), sep='_'), paste('fba', rownames(a), sep='_'), paste('v_bounds', colnames(a), sep='_'), 'objs.eql', 'gdls')
    
    #cols
    obj <- lb <- numeric(ncol(Amat)); ub <- rep(Inf, ncol(Amat))
    vtype <- rep('C', ncol(Amat))
    names(obj) <- names(lb) <- names(ub) <- names(vtype) <- colnames(Amat)
    #beta<=1
    ub[2*nrxns+nmets + 1:nrxns] <- 1
    vtype[2*nrxns+nmets + 1:nrxns] <- 'B'
    #mu free
    lb[nrxns + 1:(nmets+nrxns)] <- -Inf
    
    #rows
    rhs <- numeric(nrow(Amat)); sense <- rep('E', nrow(Amat)); names(rhs) <- names(sense) <- rownames(Amat)
    rhs['dual_eq_biomass'] <- -1; rhs[nrxns + 1:nrxns] <- lam.ub; rhs['gdls'] <- nrxns-n.ko
    sense[nrxns + 1:nrxns] <- 'L'; sense[2*nrxns+nmets + 1:nrxns] <- 'L'; sense['gdls'] <- 'G'
    
    dum <- matrix(NA, nrow=4, ncol=length(rxns.ss), dimnames=list(c('milp.obj', 'fba.obj', 'n.ko', 'stat'), rxns.ss))
    dum.lst <- list()
    for (rxn.tmp in rxns.ss){
        cat(rxn.tmp, '\n')
        rhs[paste('v_bounds', rxn.tmp, sep='_')] <- lb[rxn.tmp] <- obj[paste('l', rxn.tmp, sep='_')] <- Amat['objs.eql', paste('l', rxn.tmp, sep='_')] <- 1
        ub[paste('beta', rxn.tmp, sep='_')] <- 0
        #loose bounds on this dual
        rhs[paste('l_bounds', rxn.tmp, sep='_')] <- 10**5
        
        rc <- Rcplex(Amat=Amat, sense=sense, bvec=rhs, cvec=obj, lb=lb, ub=ub, vtype=vtype, objsense='max', control=ctrl)
        cat('milp status:', rc$status, '; objective:', rc$obj, '\n')
        v <- rc$x[1:nrxns]; mu <- rc$x[nrxns + 1:nmets]; l <- rc$x[nrxns+nmets + 1:nrxns]; betas <- rc$xopt[2*nrxns+nmets + 1:nrxns]
        names(v) <- names(l) <- names(betas) <- colnames(a); names(mu) <- rownames(a)
        
        #if feasible
        dum['stat', rxn.tmp] <- rc$status
        if (rc$status %in% c(101:102, 107)){ 
            dum['milp.obj', rxn.tmp] <- round(rc$obj, digits=4)
            dum['n.ko', rxn.tmp] <- sum(l>10**-6 & betas<10**-6)-1
            #cheating betas are 10**-3 or 1-10**-3
            cat('n cheating betas:', sum(betas<1 & !(betas %in% c(0,1))), '; n cheating lambdas:', sum(l>lam.ub*(1-betas)+10**-6), '\n')
            #fba
            fu <- rep(Inf, ncol(a)); fl <- numeric(ncol(a)); names(fu) <- names(fl) <- colnames(a)
            fu[rxn.tmp] <- fl[rxn.tmp] <- 1
            kos <- setdiff(names(betas)[l>10**-6 & betas<10**-6], rxn.tmp)
            f.tmp <- FBA(a, fba.ub=fu, fba.lb=fl, ko=kos, eps=0, control=list(trace=0))
            if (f.tmp$status==1) dum['fba.obj', rxn.tmp] <- round(f.tmp$xopt['biomass'], digits=4)
            
            #list struct for dual rxn sets
            if (f.tmp$obj>10**-6 & f.tmp$obj<median(lam.ub)){ dum.lst[[rxn.tmp]] <- setdiff(names(betas)[l>10**-6 & betas<10**-6], rxn.tmp) }
        }
        cat('\n')
        
        #reset
        rhs[paste('v_bounds', rxn.tmp, sep='_')] <- lb[rxn.tmp] <- obj[paste('l', rxn.tmp, sep='_')] <- Amat['objs.eql', paste('l', rxn.tmp, sep='_')] <- 0
        ub[paste('beta', rxn.tmp, sep='_')] <- 1
    }
    return(list(mat=dum, lst=dum.lst))
}


#dual.ub.milp w/ cut.met.mat
dum.cut <- function(a, rxns.ss=colnames(a), n.ko=Inf, ctrl=list(tilim=5), fba.ub=rep(10**5, ncol(a)), lam.ub=rep(10**3, ncol(a)), cut.met.mat=Matrix(0, nrow=1, ncol=ncol(a))){

    nrxns <- ncol(a); nmets <- nrow(a)
    fba.obj <- Matrix(0, nrow=1, ncol=nrxns, dimnames=list(1, colnames(a))); fba.obj[1,'biomass'] <- 1
    names(fba.ub) <- names(lam.ub) <- colnames(a)
    trans.rxns <- grep('^TRANS-RXN', colnames(a), value=TRUE); fba.ub[trans.rxns] <- 10**3 * mean(fba.ub)
    
    #decision vars are: v mu lambda beta
    #row blocks are: s'mu-lambda=-c, lambda+M*beta<=M, sv=0, v-M*beta<=0, dual.obj=primal.obj, cut.met.mat, gdls
    Amat <- rBind(cBind(Diagonal(x=numeric(nrxns)), t(a), -1*Diagonal(n=nrxns), Diagonal(x=numeric(nrxns))),
    cBind(Matrix(0, nrow=nrxns, ncol=nrxns+nmets), Diagonal(n=nrxns), Diagonal(x=lam.ub)),
    cBind(a, Matrix(0, nrow=nmets, ncol=2*nrxns+nmets)),
    cBind(Diagonal(n=nrxns), Matrix(0, nrow=nrxns, ncol=nmets+nrxns), Diagonal(x=-fba.ub)),
    cBind(-fba.obj, Matrix(0, nrow=1, ncol=2*nrxns+nmets)),
    cBind(Matrix(0, nrow=nrow(cut.met.mat), ncol=nrxns+nmets+nrxns), cut.met.mat),
    cBind(Matrix(rep(c(0, 1), times=c(2*nrxns+nmets, nrxns)), nrow=1, ncol=3*nrxns+nmets)))
    
    #names
    colnames(Amat) <- c(colnames(a), paste('mu', rownames(a), sep='_'), paste('l', colnames(a), sep='_'), paste('beta', colnames(a), sep='_'))
    rownames(Amat) <- c(paste('dual_eq', colnames(a), sep='_'), paste('l_bounds', colnames(a), sep='_'), paste('fba', rownames(a), sep='_'), paste('v_bounds', colnames(a), sep='_'), 'objs.eql', 
    paste('cut_met', 1:nrow(cut.met.mat), sep=''), 'gdls')
    
    #cols
    obj <- lb <- numeric(ncol(Amat)); ub <- rep(Inf, ncol(Amat))
    vtype <- rep('C', ncol(Amat))
    names(obj) <- names(lb) <- names(ub) <- names(vtype) <- colnames(Amat)
    #beta<=1
    ub[2*nrxns+nmets + 1:nrxns] <- 1
    vtype[2*nrxns+nmets + 1:nrxns] <- 'B'
    #mu free
    lb[nrxns + 1:(nmets+nrxns)] <- -Inf
    
    #rows
    rhs <- numeric(nrow(Amat)); sense <- rep('E', nrow(Amat)); names(rhs) <- names(sense) <- rownames(Amat)
    rhs['dual_eq_biomass'] <- -1; rhs[paste('l_bounds', colnames(a), sep='_')] <- lam.ub; rhs['gdls'] <- nrxns-n.ko
    sense[c(paste('l_bounds', colnames(a), sep='_'), paste('v_bounds', colnames(a), sep='_'), paste('cut_met', 1:nrow(cut.met.mat), sep=''))] <- 'L'; sense['gdls'] <- 'G'
    
    dum <- matrix(NA, nrow=4, ncol=length(rxns.ss), dimnames=list(c('milp.obj', 'fba.obj', 'n.ko', 'stat'), rxns.ss))
    dum.lst <- list()
    for (rxn.tmp in rxns.ss){
        cat(rxn.tmp, '\n')
        #set obj elements and flux/beta bounds of rxn.tmp to 1
        obj[paste('l', rxn.tmp, sep='_')] <- Amat['objs.eql', paste('l', rxn.tmp, sep='_')] <- 1
        lb[paste('beta', rxn.tmp, sep='_')] <- lb[rxn.tmp] <- ub[rxn.tmp] <- 1
        #loose bounds on this dual
        Amat[paste('l_bounds', rxn.tmp, sep='_'), paste('beta', rxn.tmp, sep='_')] <- 0
        
        rc <- Rcplex(Amat=Amat, sense=sense, bvec=rhs, cvec=obj, lb=lb, ub=ub, vtype=vtype, objsense='max', control=ctrl)
        cat('milp status:', rc$status, '; objective:', rc$obj, '\n')
        v <- rc$x[1:nrxns]; mu <- rc$x[nrxns + 1:nmets]; l <- rc$x[nrxns+nmets + 1:nrxns]; betas <- rc$xopt[2*nrxns+nmets + 1:nrxns]
        names(v) <- names(l) <- names(betas) <- colnames(a); names(mu) <- rownames(a)
        
        #if feasible
        dum['stat', rxn.tmp] <- rc$status
        if (rc$status %in% c(101:102, 107)){ 
            dum['milp.obj', rxn.tmp] <- round(rc$obj, digits=4)
            dum['n.ko', rxn.tmp] <- nrxns-sum(betas>10**-6)
            #cheating betas are 10**-3 or 1-10**-3
            cat('n cheating betas:', sum(betas<1 & !(betas %in% c(0,1))), '; n cheating lambdas:', sum(l>lam.ub*(1-betas)+10**-6), '\n')
            #fba
            fu <- rep(Inf, ncol(a)); fl <- numeric(ncol(a)); names(fu) <- names(fl) <- colnames(a)
            fu[rxn.tmp] <- fl[rxn.tmp] <- 1
            kos <- setdiff(names(betas)[betas<10**-6], rxn.tmp)
            f.tmp <- FBA(a, fba.ub=fu, fba.lb=fl, ko=kos, eps=0, control=list(trace=0))
            if (f.tmp$status==1) dum['fba.obj', rxn.tmp] <- round(f.tmp$xopt['biomass'], digits=4)
            
            #list struct for dual rxn sets
            if (f.tmp$obj>10**-6 & f.tmp$obj<median(lam.ub)){ dum.lst[[rxn.tmp]] <- setdiff(names(l)[l>10**-6 & betas<10**-6], rxn.tmp) }
        }
        cat('\n')
        
        #reset
        obj[paste('l', rxn.tmp, sep='_')] <- Amat['objs.eql', paste('l', rxn.tmp, sep='_')] <- 0
        lb[paste('beta', rxn.tmp, sep='_')] <- lb[rxn.tmp] <- 0
        ub[rxn.tmp] <- Amat[paste('l_bounds', rxn.tmp, sep='_'), paste('beta', rxn.tmp, sep='_')] <- fba.ub[rxn.tmp]
        
    }
    return(list(mat=dum, lst=dum.lst))
}
