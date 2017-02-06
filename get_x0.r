##jmd
##2.22.12
##get_x0.r
##get beta that optimizes probabilities and given growth conditions

#decision vars are [v1 v2 ... accum_v deplete_v accum.slack deplete.slack beta], where v1 & v2 are for metab.dilution constraints
#amat=[S 0 0 0 0; 0 S_intra 0 -delta*|S_intra|; 0 0 S_intra delta*|S_intra|; I 0 0 -diag(ub); 0 I 0 -diag(md.ub); 0 0 I -diag(md.ub)]
#'goals' are treated as requirements here, inconsistent w/ FBA()
#f.ub.mat has rows per flux and columns per condition
#bm.ub st if bm_flux is a rewarded goal, it won't try to stretch unrealistically high
get.x0 <- function(sp, obj.coeffs, goal.coeff=10**3, check.input=TRUE, se='E', ctrl=list(tilim=120), coeff.mets=10**-3, biomass='biomass', 
nondilute.regexp='^(PROTON|WATER|OH|OXYGEN-MOLECULE|CARBON-DIOXIDE)(\\[CCO-|$)|\\[CCO-EXTRACELLULAR\\]$', known.rxns=NULL, ko.rxns.lst=NULL,
md=TRUE, connect.constr.mat=NULL, mipstart=NULL, 
f.ub.mat=array(10**3, dim=ncol(sp), dimnames=list(colnames(sp))), 
md.f.ub=array(10**3, dim=ncol(sp), dimnames=list(colnames(sp))),
bm.ub=0.4, req.conds=ncol(as.matrix(f.ub.mat)), lb.req=0.1, goal.conds=setdiff(ncol(f.ub.mat), req.conds) ){
    
    f.ub.mat <- as.matrix(f.ub.mat)
    
    #make 3col amat
    sp.3c <- dense2sp(sp)
    
    #get scalars
    nconds <- ncol(f.ub.mat)
    nrxns <- ncol(sp)
    nmets <- nrow(sp)
    
    se.v <- rep(NA, nmets)
    se.v[1:nmets] <- se
        
    if (check.input){
        stopifnot(!is.na(obj.coeffs), all(!is.na(sp)), !is.na(f.ub.mat), is.finite(f.ub.mat), nrow(as.matrix(f.ub.mat))==ncol(sp),
        is.null(names(se))|all(names(se)==rownames(sp)), is.null(names(obj.coeffs))|all(names(obj.coeffs)==colnames(sp)), 
        is.null(names(coeff.mets))|all(names(coeff.mets)==rownames(sp)), is.null(connect.constr.mat)|sum(duplicated(rownames(connect.constr.mat)))==0
        )
        cat('Done input validity checks.\n')
    }
    
    amat.3c <- NULL
    
    ##for each condition, rbind an S matrix for Sv=0, & an identity matrix and -Diagonal(f.ub.mat[,i]) matrix for v<=f.ub*beta
    for (cond.ind in 1:nconds){
        amat.3c <- rbind(amat.3c, 
        data.frame(rxn=paste('cond', cond.ind, '_', sp.3c$rxn, sep=''), cpd=paste('cond', cond.ind, '_', sp.3c$cpd, sep=''), coeff=sp.3c$coeff),
        data.frame(rxn=paste('cond', cond.ind, '_', colnames(sp), sep=''), cpd=paste('cond', cond.ind, '_', colnames(sp), sep=''), coeff=1),
        data.frame(rxn=paste('beta_', colnames(sp), sep=''), cpd=paste('cond', cond.ind, '_', colnames(sp), sep=''), coeff=-f.ub.mat[,cond.ind])
        )
    }
    #name cols and rows of Sv=0 and v<=ub*beta (these rows have same name as cols)
    cnames <- paste(rep(paste('cond', 1:nconds, '_', sep=''), each=nrxns), rep(colnames(sp), times=nconds), sep='')
    obj <- c(numeric(nconds*nrxns))
    ub <- as.numeric(f.ub.mat)
    lb <- numeric(length(ub))
    
    rnames <- c(paste(rep(paste('cond', 1:nconds, '_', sep=''), each=nmets), rep(rownames(sp), times=nconds), sep=''), cnames)
    sense <- c(rep(se.v, times=nconds), rep('L', nrxns*nconds))
    rhs <- numeric(length(rnames))
    
    ##md constraints
    if (md){
        mets.intra <- grep(nondilute.regexp, rownames(sp)); nmets.intra <- length(mets.intra)
        
        md.res <- run.get.md.constraints.3c(sp=sp, mets.intra=mets.intra, md.f.ub=md.f.ub)
        amat.3c <- rbind(amat.3c, md.res$amat)
        cnames <- c(cnames, md.res$cnames)
        rnames <- c(rnames, md.res$rnames)

        obj <- c(obj, numeric(2*nrxns), rep(-10**3, 2*nmets.intra))
        ub <- c(ub, md.f.ub, md.f.ub, rep(1000, 2*nmets.intra))
        lb <- numeric(length(ub))
        rhs <- numeric(length(rnames))
        sense <- c(sense, rep('L', length(md.res$rnames)))
    }
        
    ##dual.constraints
    if (!is.null(ko.rxns.lst) & length(ko.rxns.lst)>0){
        dual.res <- get.dual.constraints(sp=sp, ko.rxns.lst=ko.rxns.lst, biomass=biomass, keep.rxns=known.rxns)
        amat.3c <- rbind(amat.3c, dual.res$amat.3c)
        
        cnames <- c(cnames, names(dual.res$ub))
        lb <- c(lb, dual.res$lb)
        ub <- c(ub, dual.res$ub)
        obj <- c(obj, dual.res$obj)
        
        rhs <- c(rhs, dual.res$rhs)
        sense <- c(sense, dual.res$sense)
        rnames <- c(rnames, names(dual.res$rhs))
    }
        
    ##connect.constr.mat
    if (!is.null(connect.constr.mat)){
        ccm.3c <- dense2sp(connect.constr.mat)
        amat.3c <- rbind(amat.3c, data.frame(rxn=paste('beta_', ccm.3c$rxn, sep=''), cpd=paste('ccm_', ccm.3c$cpd, sep=''), coeff=ccm.3c$coeff))
        rnames <- c(rnames, paste('ccm_', rownames(connect.constr.mat), sep=''))
        sense <- c(sense, rep('L', nrow(connect.constr.mat)))
        rhs <- c(rhs, numeric(nrow(connect.constr.mat)))
    }
    
    ##add beta: obj=obj.coeffs, ub=1
    obj <- c(obj, obj.coeffs)
    ub <- c(ub, rep(1, nrxns))
    lb <- c(lb, numeric(nrxns))
    cnames <- c(cnames, paste('beta_', colnames(sp), sep=''))
    names(obj) <- names(lb) <- names(ub) <- cnames
    
    ##names
    rownames.f <- factor(amat.3c$cpd, levels=rnames, ordered=TRUE)
    colnames.f <- factor(amat.3c$rxn, levels=cnames, ordered=TRUE)
    
    ##amat
    amat <- sparseMatrix(i=as.numeric(rownames.f), j=as.numeric(colnames.f), x=amat.3c$coeff, dimnames=list(levels(rownames.f), levels(colnames.f)))
    
    ##vector names
    stopifnot(rownames(amat)==rnames, colnames(amat)==cnames)
    names(sense) <- names(rhs) <- rownames(amat)
    
    ##goals/reqs
    if (length(goal.conds)>0){ 
        obj[paste('cond', goal.conds, '_', biomass, sep='')] <- goal.coeff 
        ub[paste('cond', goal.conds, '_', biomass, sep='')] <- bm.ub
    }
    if (!is.null(req.conds)){ lb[paste('cond', req.conds, '_', biomass, sep='')] <- lb.req }
    if (!is.null(known.rxns)){ lb[paste('beta_', known.rxns, sep='')] <- 1 }
    
    ##rcplex
    if (is.null(mipstart)){ vtype <- rep('C', ncol(amat)) } else { vtype <- rep(c('C', 'B'), times=c(ncol(amat)-nrxns, nrxns)) }
    #(mi)lp
    mp <- Rcplex(cvec=obj, Amat=amat, bvec=rhs, ub=ub, lb=lb, sense=sense, objsense='max', vtype=vtype, mipstart=mipstart, control=ctrl)
    cat('CPLEX status of LP', mp$stat, '\n')
    #return rcplex
    x <- mp$xopt
    names(x) <- names(mp$xopt) <- names(obj)
    beta.names <- grep('^beta_', names(x), value=TRUE)
    #names(x0) <- gsub('^beta_', '', x=names(x0))
    accum.slack <- x[grep('^accum_slack_', names(x))]
    deplete.slack <- x[grep('^deplete_slack_', names(x))]
    
    ##what to return? beta.v & accum/deplete slack
    ret.mat <- matrix(cbind(x[beta.names], obj[beta.names]), ncol=2, dimnames=list(gsub('^beta_', '', beta.names), c('beta.v', 'obj.coeff')))
    ret.rcplex <- mp
    if (!md){ ret.md <- NULL } else { ret.md <- matrix(cbind(accum.slack, deplete.slack), ncol=2, dimnames=list(mets.intra, c('accum', 'deplete'))) }
    return(list(mat=ret.mat, md=ret.md, rcplex=ret.rcplex))
}
