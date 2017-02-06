##jmd
##2.27.12
##get_dual_constraints.r

#amat has columns for beta, but ub/lb/obj do not have elements for beta

get.dual.constraints <- function(sp, ko.rxns.lst, biomass='biomass', keep.rxns=NULL){
    n.dual.conds <- length(ko.rxns.lst)
    nrxns <- ncol(sp)
    nmets <- nrow(sp)
    
    dual.ub <- rep(1, nrxns)
    names(dual.ub) <- colnames(sp)
    keep.rxns <- setdiff(keep.rxns, biomass)
    dual.ub[intersect(keep.rxns, names(dual.ub))] <- 0
    
    ub <- rep(Inf, (nrxns+nmets)*n.dual.conds)
    lb <- rep(-Inf, (nrxns+nmets)*n.dual.conds)
    obj <- numeric((nrxns+nmets)*n.dual.conds)
    names(ub) <- names(lb) <- names(obj) <- paste('dual_cond', rep(1:n.dual.conds, each=nrxns+nmets), '_', rep(c(rownames(sp), colnames(sp)), times=n.dual.conds), sep='')
    
    fba.obj <- named.vec(0, colnames(sp)); fba.obj[biomass] <- 1
    
    sense <- rep(rep(c('E', 'L'), each=nrxns), times=n.dual.conds)
    rhs <- rep(c(-fba.obj, dual.ub), times=n.dual.conds)
    names(rhs) <- names(sense) <- paste('dual_cond', rep(1:n.dual.conds, each=2*nrxns), '_', rep(c('eq', 'ko'), each=nrxns),  '_', colnames(sp), sep='')
                    
    ##bounds & penalties on nut rxns
    bm.cols <- paste('dual_cond', 1:n.dual.conds, '_', biomass, sep='')
    stopifnot(bm.cols %in% names(lb))
    lb[bm.cols] <- 0
    obj[bm.cols] <- 100
    
    amat.3c <- NULL
    #decision vars are: [m1, r1, m2, r2, ..., beta]
    for (dual.cond.i in 1:n.dual.conds){
        ko.rxns <- ko.rxns.lst[[dual.cond.i]]
                
        ##assign betas
        #dci.beta.coeff.v is coeff of beta in these rows of amat
        dci.beta.coeff.v <- dual.ub
        #don't want beta=1 -> lambda<0 for ko.rxns or nutrient rxns
        rxns.no.beta <- c(biomass, ko.rxns)
        dci.beta.coeff.v[intersect(names(dci.beta.coeff.v), rxns.no.beta)] <- 0
        
        #decision vars for this case are [m_i, r_i, beta]
        amat.tmp <- rBind(cBind(t(sp), -1*Diagonal(n=nrxns), Matrix(0, nrow=nrxns, ncol=nrxns)),
        cBind(Matrix(0, nrow=nrxns, ncol=nmets), Diagonal(n=nrxns), Diagonal(x=dci.beta.coeff.v)))
        
        rownames(amat.tmp) <- paste('dual_cond', dual.cond.i, '_', rep(c('eq', 'ko'), each=nrxns),  '_', colnames(sp), sep='')
        colnames(amat.tmp) <- c(paste('dual_cond', dual.cond.i, '_', c(rownames(sp), colnames(sp)), sep=''), paste('beta_', colnames(sp), sep=''))
                
        amat.tmp.3c <- dense2sp(amat.tmp)
        amat.3c <- rbind(amat.3c, amat.tmp.3c)
    }#end for
    return(list(amat.3c=amat.3c, obj=obj, ub=ub, lb=lb, sense=sense, rhs=rhs))
}
