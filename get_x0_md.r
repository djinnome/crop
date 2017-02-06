##jmd
##2.7.12
##get_x0_multi.r

#for single growth scenario w/ md constraints
#decision vars are [v1 v2 v3 accum.slack deplete.slack beta], where v1 & v2 are for metab.dilution constraints
#amat=[S 0 0 0 0; 0 S_intra 0 -delta*|S_intra|; 0 0 S_intra delta*|S_intra|; I 0 0 -diag(ub); 0 I 0 -diag(md.ub); 0 0 I -diag(md.ub)]
#'goals' are treated as requirements here, inconsistent w/ FBA()
get.x0.md <- function(sp, obj.coeffs, goals='biomass', lb.goals=1, check.input=TRUE, se='E', mipstart=NULL, ctrl=list(tilim=120), coeff.mets=10**-3,
nondilute.regexp='^(PROTON|WATER|OH|OXYGEN-MOLECULE|CARBON-DIOXIDE)(\\[CCO-|$)|\\[CCO-EXTRACELLULAR\\]$', 
connect.constr.mat=matrix(nrow=0, ncol=ncol(sp), dimnames=list(NULL, colnames(sp))), 
f.ub=array(10**3, dim=ncol(sp), dimnames=list(colnames(sp))), 
md.f.ub=array(10**3, dim=ncol(sp), dimnames=list(colnames(sp))) ){
    
    #make sp.intra, which has same number of rxns as sp
    sp.intra <- sp[-grep(nondilute.regexp, rownames(sp)),]
    sp.intra.norm <- abs(sp.intra)/rowSums(abs(sp.intra))
    
    nrxns <- ncol(sp)
    nmets <- nrow(sp)
    nmets.intra <- nrow(sp.intra)
        
    if (check.input){
        stopifnot(
        goals %in% colnames(sp), !is.na(obj.coeffs), all(!is.na(sp)), !is.na(f.ub), is.finite(f.ub), !is.null(names(f.ub)), ncol(sp)==ncol(sp.intra),
        is.null(mipstart)|length(mipstart)==4*ncol(sp)+2*nmets.intra,  
        is.null(names(se))|all(names(se)==rownames(sp)), is.null(names(obj.coeffs))|all(names(obj.coeffs)==colnames(sp)), 
        is.null(names(coeff.mets))|all(names(coeff.mets)==rownames(sp)) 
        )
        cat('Done input validity checks.\n')
    }
    
#   if(length(s.met.coeffs)==1){ 
#        s.met.coeffs.v <- array(data=s.met.coeffs, dim=nrow(sp.intra), dimnames=rownames(sp.intra)) 
#    } else { 
#        s.met.coeffs.v <- array(data=s.met.coeffs[rownames(sp.intra)], dim=nrow(sp.intra), dimnames=rownames(sp.intra))
#    }#end else
    
    amat <- rBind(
    cBind(sp, Matrix(0, nr=nmets, nc=3*nrxns+2*nmets.intra)),
    cBind(Matrix(0, nr=nmets.intra, nc=nrxns), -1*sp.intra, Matrix(0, nr=nmets.intra, nc=nrxns), -1*Diagonal(n=nmets.intra), Matrix(0, nr=nmets.intra, nc=nmets.intra), coeff.mets * sp.intra.norm),
    cBind(Matrix(0, nr=nmets.intra, nc=2*nrxns), sp.intra, Matrix(0, nr=nmets.intra, nc=nmets.intra), -1*Diagonal(n=nmets.intra), coeff.mets * sp.intra.norm),
    cBind(Diagonal(n=nrxns), Matrix(0, nr=nrxns, nc=2*nrxns+2*nmets.intra), Diagonal(x=-f.ub)),
    cBind(Matrix(0, nr=nrxns, nc=nrxns), Diagonal(n=nrxns), Matrix(0, nr=nrxns, nc=nrxns+2*nmets.intra), Diagonal(x=-md.f.ub)),
    cBind(Matrix(0, nr=nrxns, nc=2*nrxns), Diagonal(n=nrxns), Matrix(0, nr=nrxns, nc=2*nmets.intra), Diagonal(x=-md.f.ub)),
    cBind(Matrix(0, nr=nrow(connect.constr.mat), nc=3*nrxns+2*nmets.intra), connect.constr.mat)
    )

    obj <- c(numeric(3*nrxns), rep(-10**4, 2*nmets.intra), obj.coeffs)
    sense <- c(array(data=se, dim=nrow(sp)), rep('L', 2*nmets.intra+3*nrxns+nrow(connect.constr.mat)))
    rhs <- numeric(nrow(amat))
    lb <- numeric(ncol(amat))
    ub <- c(f.ub, md.f.ub, md.f.ub, rep(1000, 2*nmets.intra), rep(1, nrxns))
    
    rownames(amat) <- names(sense) <- names(rhs) <- c(
    paste(rep(c('', 'accum_', 'deplete_'), times=c(nmets, nmets.intra, nmets.intra)), c(rownames(sp), rownames(sp.intra), rownames(sp.intra)), sep=''), 
    paste(rep(c('', 'accum_', 'deplete_'), each=nrxns), rep(colnames(sp), times=3), sep=''), 
    rownames(connect.constr.mat)) 
    
    colnames(amat) <- names(lb) <- names(ub) <- names(obj) <- paste(rep(c('', 'accum_', 'deplete_', 'accum_slack_', 'deplete_slack_', 'beta_'), times=c(nrxns, nrxns, nrxns, nmets.intra, nmets.intra, nrxns)),
    c(rep(colnames(sp), times=3), rownames(sp.intra), rownames(sp.intra), colnames(sp)), sep='')
    
    ##if infeasible
    #don't really need depletion, since can add exports later
    #rhs[paste('deplete_', rownames(sp.intra), sep='')] <- 10**3
    #rhs[paste('accum_', rownames(sp.intra), sep='')] <- -Inf
    #try w/ ub[3*nrxns+1:nrxns]=lb[3*nrxns+1:nrxns]=as.numeric(colnames(sp.f) %in% rxns1.unblocked), so can see inconsistent rxns
    
    #goals lb
    lb[goals] <- lb.goals
    lb[paste('beta_', goals, sep='')] <- 1
    if (is.null(mipstart)){ vtype <- rep('C', ncol(amat)) } else { vtype <- rep(c('C', 'B'), times=c(ncol(amat)-nrxns, nrxns)) }
    #(mi)lp
    mp <- Rcplex(cvec=obj, Amat=amat, bvec=rhs, ub=ub, lb=lb, sense=sense, objsense='max', vtype=vtype, mipstart=mipstart, control=ctrl)
    cat('CPLEX status of LP', mp$stat, '\n')
 
    x <- mp$xopt
    names(x) <- names(obj)
    v <- x[1:nrxns]
    x0 <- x[grep('^beta_', names(x))]
    #names(x0) <- gsub('^beta_', '', x=names(x0))
    accum.slack <- x[grep('^accum_slack_', names(x))]
    deplete.slack <- x[grep('^deplete_slack_', names(x))]
    
    return(
    list(mat=matrix(cbind(v, x0, obj[1:nrxns]), ncol=3, dimnames=list(names(obj)[1:nrxns], c('v', 'beta.v', 'obj.coeff'))), 
    slack=matrix(cbind(accum.slack, deplete.slack), ncol=2, dimnames=list(rownames(sp.intra), c('accum', 'deplete'))))
    )
}
