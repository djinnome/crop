##jmd
##2.7.12
##get_x0_simple.r

#for single growth scenario w/o md constraints
#decision vars are [v beta]
#obj.sense=MAX
#'goals' are treated as requirements
get.x0.simple <- function(sp, obj.coeffs, goals='biomass', lb.goals=1, check.input=TRUE, se='E', mipstart=NULL, ctrl=list(tilim=60),
connect.constr.mat=matrix(nrow=0, ncol=ncol(sp), dimnames=list(NULL, colnames(sp))), f.ub=array(10**3, dim=ncol(sp), dimnames=list(colnames(sp))) ){
    if (check.input){
        stopifnot(names(obj.coeffs)==colnames(sp), goals %in% colnames(sp), !is.na(obj.coeffs), all(!is.na(sp)), !is.na(f.ub), is.finite(f.ub), !is.null(names(f.ub)), 
        is.null(names(se))|all(names(se)==rownames(sp)), is.null(mipstart)|length(mipstart)==2*ncol(sp))
        cat('Done input validity checks.\n')
    }
    
    obj <- c(numeric(ncol(sp)), obj.coeffs)
    fba.se <- array(data=se, dim=nrow(sp), dimnames=rownames(sp))
    #A=[s 0; I -f.ub]; decision vars are [v beta]
    amat <- rBind(cBind(sp, Matrix(0, nrow=nrow(sp), ncol=ncol(sp))),
    cBind(Diagonal(n=ncol(sp)), -1*Diagonal(x=f.ub)))
    
    #add connect.constr.mat
    amat <- rBind(amat, cBind(Matrix(0, nrow=nrow(connect.constr.mat), ncol=ncol(sp)), connect.constr.mat))
    sense <- c(fba.se, rep('L', times=ncol(sp)), rep('L', nrow(connect.constr.mat)))
    rownames(amat) <- c(rownames(sp), colnames(sp), rownames(connect.constr.mat))

    rhs <- numeric(nrow(amat))
    names(sense) <- names(rhs) <- rownames(amat)
    
    ub <- c(f.ub, rep(1, ncol(sp)))
    lb <- numeric(ncol(amat))
    colnames(amat) <- names(lb) <- names(ub) <- names(obj) <- c(colnames(sp), paste('beta_', colnames(sp), sep=''))
    #goals lb
    lb[goals] <- lb.goals
    lb[paste('beta_', goals, sep='')] <- 1
    if (is.null(mipstart)){ vtype <- rep('C', ncol(amat)) } else { vtype <- rep(c('C', 'B'), each=ncol(sp)) }
    #(mi)lp
    mp <- Rcplex(cvec=obj, Amat=amat, bvec=rhs, ub=ub, lb=lb, sense=sense, objsense='max', vtype=vtype, mipstart=mipstart, control=ctrl)
    cat('CPLEX status of LP', mp$stat, '\n')
    x <- mp$xopt
    names(x) <- names(obj)
    v <- x[1:ncol(sp)]
    x0 <- x[ncol(sp) + 1:ncol(sp)]
    names(x0) <- gsub('beta_', '', x=names(x0))
    
    return(list(v=v, beta.v=x0))
}
