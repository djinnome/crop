##jmd
##8.22.11
##min_met_adj.r

#checks for input/output of mets
#penalize is set of reactions and/or metabolites that get very high penalty
min.met.adj <- function(sp, penalize=NULL, fba.ub=rep(10**6, ncol(sp)), bm.min=1, ko.rxns=NULL){
    stopifnot(penalize %in% c(colnames(sp), rownames(sp)))
    nmets <- nrow(sp); nrxns <- ncol(sp)
    #decision vars are [v r+ r-]
    obj <- c(rep(0, nrxns), rep(c(1, 0.1), each=nmets))
    a <- cBind(sp, Diagonal(n=nmets), -1*Diagonal(n=nmets))
    lb <- rep(0, nrxns+2*nmets)
    ub <- c(fba.ub, rep(10**6, 2*nmets))

    names(obj) <- names(lb) <- names(ub) <- c(colnames(sp), rownames(sp), paste('neg', rownames(sp), sep='_'))
    if(!is.null(ko.rxns)){ ub[ko.rxns] <- 0 }
    lb['biomass'] <- bm.min
    #penalize
    pen <- c(as.character(unique(bm$comp[bm$coeff>0])), penalize)
    obj[pen] <- 10**6
    #run lp
    mfa <- Rcplex(cvec=obj, Amat=a, bvec=rep(0, nrow(a)), lb=lb, ub=ub, sense='E', objsense='min')
    names(mfa$xopt) <- names(obj)
    r <- mfa$xopt[-(1:nrxns)]
    ret <- r[r>0]
    #verify soln
    cat('Verifying soln \n')
    fba.se <- rep('E', nrow(sp)); names(fba.se) <- rownames(sp)
    fba.se[rownames(sp)[r[1:nmets]>10**-6]] <- 'L'; fba.se[rownames(sp)[r[nmets + 1:nmets]>10**-6]] <- 'G'
    f <- FBA(a=sp, control=list(trace=0), fba.ub=fva.ub, sense=fba.se)
   
    if (f$obj>10**-6) return(ret[order(-ret)]) else return(NA)
}

#checks for input/output of mets
#can weight metabolites by mass, so find smallest one
#can have n.max.import & n.max.export; can set to only export by n.max.import=0
#can have fva.ub instead of M
mma.bin <- function(sp, fva.ub=rep(1000, ncol(sp)), penalize=NULL, mets0=NULL, n.add=5, import.wts=rep(1, nrow(sp)), export.wts=import.wts, ctrl=list(trace=1, tilim=30), eps=10**-6, M=10**3, bm.lb=0.1){
    stopifnot(penalize %in% c(colnames(sp), rownames(sp)))
    nmets <- nrow(sp); nrxns <- ncol(sp)
    
    #decision vars are [v r+ r- beta+ beta-]
    #rows are [sv + r+ - r- = 0; r+ <= 10^6*beta+; r- <= 10^6*beta-; sum(new beta+)<=n.add]
    obj <- c(rep(0, nrxns+2*nmets), import.wts, export.wts)
    a <- rBind(cBind(sp, Diagonal(n=nmets), -1*Diagonal(n=nmets), Matrix(0, nrow=nmets, ncol=2*nmets)),
    cBind(Matrix(0, nrow=nmets, ncol=nrxns), Diagonal(n=nmets), Matrix(0, nrow=nmets, ncol=nmets), -M*Diagonal(n=nmets), Matrix(0, nrow=nmets, ncol=nmets)),
    cBind(Matrix(0, nrow=nmets, ncol=nrxns+nmets), Diagonal(n=nmets), Matrix(0, nrow=nmets, ncol=nmets), -M*Diagonal(n=nmets)),
    cBind(Matrix(0, nrow=1, ncol=nrxns+2*nmets), as.numeric(!(rownames(sp) %in% mets0)), Matrix(0, nrow=1, ncol=nmets)))
    lb <- numeric(length(obj))
    #Inf ub on fluxes
    ub <- rep(Inf, length(obj))
    ub[1:ncol(sp)] <- fva.ub
    ub[-(1:(nrxns+2*nmets))] <- 1
    names(obj) <- names(lb) <- names(ub) <- c(colnames(sp), rownames(sp), paste('neg', rownames(sp), sep='_'), paste('beta_', rownames(sp), sep=''), paste('beta_neg', rownames(sp), sep='_'))
    lb['biomass'] <- bm.lb
    #penalize
    pen <- paste('beta_', c(as.character(unique(bm$comp[bm$coeff>0])), penalize), sep='')
    obj[pen] <- 10**6
    #rhs
    sense <- rep(c('E', 'L'), times=c(nmets, 2*nmets+1))
    rhs <- c(numeric(nrow(a)-1), n.add)
    #vtype
    vt <- rep('C', ncol(a))
    vt[nrxns+2*nmets + 1:nmets] <- 'B'
    #run milp
    mfa <- Rcplex(cvec=obj, Amat=a, bvec=rhs, lb=lb, ub=ub, sense=sense, objsense='min', control=ctrl, vtype=vt)
    cat('MILP status', mfa$status, 'w/ obj', mfa$obj, '\n')
    names(mfa$xopt) <- names(obj)
    beta.v <- mfa$xopt[nrxns+2*nmets + 1:(2*nmets)]
    ret <- beta.v[beta.v>eps]
    ret <- ret[order(!(names(ret) %in% grep('neg_', names(ret), value=TRUE)), ret, decreasing=TRUE)]
    #verify
    cat('Verifying soln \n')
    fba.se <- rep('E', nrow(sp)); names(fba.se) <- rownames(sp)
    fba.se[rownames(sp)[beta.v[1:nmets]>eps]] <- 'L'; fba.se[rownames(sp)[beta.v[nmets + 1:nmets]>eps]] <- 'G'
    f <- FBA(a=sp, fba.ub=fva.ub, control=list(trace=0), sense=fba.se)
    
    return(ret)
}
