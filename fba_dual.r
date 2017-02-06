##jmd
##3.8.11
##fba_dual.r

##for no-growth conds
#make matrices st can still have beta as a decision vars (tho they're restricted to 1)
#ko.rxns are rxns ko'd a priori e.g. thioredoxen or o2
#mu's >= 0 for sv>=0

fba.dual <- function(a, ko.rxns=NULL, c.rxns=c('SUCROSE-TRANS-RXN-L2R','CIT-TRANS-RXN-L2R')){
    ##instantiate
    fba.obj <- rep(0, ncol(a))
    dual.lb <- rep(c(0,-Inf, 1), times=c(nrow(a), ncol(a), ncol(a)))
    dual.ub <- rep(c(Inf, 1), times=c(nrow(a)+ncol(a), ncol(a)))
    dual.obj <- rep(0, nrow(a)+ncol(a)+ncol(a))
    
    ##decision vars are [mu lambda 'beta'] of length m+n+n
    #eq: s'*mu-I*lambda=-c
    eq <- cBind(t(a), -1*Diagonal(ncol(a)), Matrix(0,ncol(a),ncol(a)))
    #lambda+1000*beta<=1000 except for sugars & ko's
    ko.mat <- cBind(Matrix(0,ncol(a),nrow(a)), Diagonal(ncol(a)), 1000*Diagonal(ncol(a)))
    ko.rhs <- rep(1000, ncol(a))
    #names
    col.names <- c(rownames(a), colnames(a), paste('beta', colnames(a), sep='_'))
    rownames(eq) <- rownames(ko.mat) <- names(ko.rhs) <- names(fba.obj) <- colnames(a)
    colnames(eq) <- colnames(ko.mat) <- names(dual.lb) <- names(dual.ub) <- names(dual.obj) <- col.names
    
    ##special cases
    fba.obj['biomass'] <- 1
    dual.obj[names(dual.obj) %in% c.rxns] <- 100
    #beta's not relevant to c.rxns or ko.rxns
    #can set coeff of lambda and/or beta to 0; set coeff of beta = 0
    ko.mat[rownames(ko.mat) %in% c(c.rxns, ko.rxns), colnames(ko.mat) %in%  paste('beta', c(c.rxns, ko.rxns), sep='_')] <- 0
    #bounds
    dual.lb[names(dual.obj) %in% c.rxns] <- 0
    dual.ub[names(dual.obj) %in% c.rxns] <- Inf
    
    ##constraint mat
    dual.amat <- rBind(eq, ko.mat)
    dual.rhs <- c(-fba.obj, ko.rhs)
    dual.sense <- rep(c('E', 'L'), c(nrow(eq), nrow(ko.mat)))
    
    ##lp
    dual.lp <- Rcplex(cvec=dual.obj, Amat=dual.amat, bvec=dual.rhs, lb=dual.lb, ub=dual.ub, sense=dual.sense)
    #print
    cat('Status:', dual.lp$status, 'w/ dual.obj:', dual.lp$obj, '\n \n')
    #examine
    x <- dual.lp$xopt; names(x) <- names(dual.obj); x[c.rxns]
    #return
    return(dual.lp)
}#end fcn
