##jmd
##2.15.11
##min_export.r
##max weights*r: Sv-r=0, r>=0

#w applies to mets
#penalize.rxns & low.pen apply only to columns of sp
min.export <- function(sp, w=rep(1, nrow(sp)), penalize.rxns=NULL, goals='biomass', v.min=1, low.pen=NULL, se.fba=rep('E', nrow(sp)), ub.fba=rep(10**6, ncol(sp)), ctrl=list(tilim=30), milp.only=TRUE){
    nrxns <- ncol(sp); nmets <- nrow(sp)
    trans.rxns <- grep('^TRANS-RXN', colnames(sp), value=TRUE)
    names(w) <- rownames(sp)
    w[penalize.rxns] <- 10**6; w[low.pen] <- 10**-3
    
    ##lp
    #decision vars are [v r]
    if (!milp.only){
        obj <- c(rep(0, nrxns), w)
        a <- cBind(sp, -1*Diagonal(n=nmets))
        colnames(a) <- c(colnames(sp), rownames(sp)); rownames(a) <- rownames(sp)
        lb <- rep(0, nrxns+nmets)
        ub <- c(ub.fba, rep(10**6, times=nmets))
        names(obj) <- names(lb) <- names(ub) <- c(colnames(sp), rownames(sp))
        lb[goals] <- v.min
        #ub[trans.rxns] <- 10**2 * median(ub)
        
        min.ex <- Rcplex(cvec=obj, Amat=a, bvec=numeric(nrow(a)), lb=lb, ub=ub, sense='E', objsense='min', control=list(trace=0))
        names(min.ex$xopt) <- names(obj)
        cat('CPLEX status of LP:', min.ex$status, '\n')
        r <- min.ex$xopt[-(1:nrxns)]
    }
        
    ##milp
    #decision vars are [v r beta], where r is the met export flux; a=[s -I 0; 0 I -1000*I]. ncol=nrxns+nmets+nmets, nrow=nmets+nmets
    obj2 <- c(rep(0, nrxns+nmets), w)
    lb2 <- numeric(nrxns+2*nmets)
    ub2 <- c(ub.fba, rep(c(10**3, 1), times=c(nmets, nmets)))
    names(obj2) <- names(lb2) <- names(ub2) <- c(colnames(sp), rownames(sp), paste('beta', rownames(sp), sep='_'))
    lb2[goals] <- v.min
    #ub2[trans.rxns] <- 10**8
    if (!milp.only){ ub2[paste('beta', rownames(sp), sep='_')][r==0] <- 0 }
    vt <- rep(c('C', 'B'), times=c(nrxns+nmets, nmets))
    sense2 <- c(se.fba, rep('L', nmets))
    #a
    a2 <- rBind(cBind(sp, -1*Diagonal(n=nmets), Matrix(0, nrow=nmets, ncol=nmets)),
    cBind(Matrix(0, nrow=nmets, ncol=nrxns), Diagonal(n=nmets), -10^3*Diagonal(n=nmets)))
    colnames(a2) <- names(obj2); rownames(a2) <- c(rownames(sp), paste('r', rownames(sp), sep='_'))
        
    min.ex.milp <- Rcplex(cvec=obj2, Amat=a2, bvec=numeric(nrow(a2)), lb=lb2, ub=ub2, sense=sense2, vtype=vt, objsense='min', control=ctrl)
    names(min.ex.milp$xopt) <- names(obj2)
    cat('CPLEX status of MILP:', min.ex.milp$status, '\n')
    r2 <- min.ex.milp$xopt[nrxns + 1:nmets]
    
    ##return
    if (min.ex.milp$status %in% c(101,102)){
        ret <- r2[r2>0]
        print(ret[order(-ret)])
        
        ##verify soln
        if (length(ret)>=1){
            cat('Verifying soln \n')
            names(se.fba) <- rownames(sp); se.fba[names(ret)] <- 'G'
            f <- FBA(a=sp, fba.ub=ub.fba, control=list(trace=0), sense=se.fba)
        }
    }#end milp status
    
    #ex <- names(r)[r>0]
    #cat('Exports:', paste(ex, collapse=', '), '\n')
    return(ret)
}#end fcn
