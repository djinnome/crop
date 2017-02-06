##jmd
##9.27.11
##min_pen.r

min.pen <- function(sp, pen.v, goals='biomass', f.ub=rep(10**3, ncol(sp)), lb.goals=1, se=rep('E', nrow(sp)), check.input=TRUE){
    if (check.input){
        stopifnot(names(pen.v)==colnames(sp), goals %in% colnames(sp), !is.na(pen.v), all(!is.na(sp)), !is.na(f.ub), is.finite(f.ub))
        cat('Done input validity checks.\n')
    }
    
    obj <- c(numeric(ncol(sp)), pen.v)
    #A=[s 0; I -f.ub]; decision vars are [v beta]
    amat <- rBind(cBind(sp, Matrix(0, nrow=nrow(sp), ncol=ncol(sp))), 
    cBind(Diagonal(n=ncol(sp)), -1*Diagonal(x=f.ub)))
    rhs <- numeric(nrow(sp)+ncol(sp))
    ub <- rep(c(Inf, 1), each=ncol(sp))
    lb <- numeric(2*ncol(sp))
    #1st sense is for Sv (>)= 0
    sense <- c(se, rep('L', times=ncol(sp)))
    #names
    rownames(amat) <- names(sense) <- names(rhs) <- c(rownames(sp), colnames(sp))
    colnames(amat) <- names(lb) <- names(ub) <- names(obj) <- c(colnames(sp), paste('beta_', colnames(sp), sep=''))
    #lb
    lb[goals] <- lb.goals
    #lp
    mp <- Rcplex(cvec=obj, Amat=amat, bvec=rhs, ub=ub, lb=lb, sense=sense)
    cat('CPLEX status of LP', mp$stat, '\n')
    x <- mp$xopt; names(x) <- names(obj)
    v <- x[1:ncol(sp)]
    beta.v <- x[ncol(sp) + 1:ncol(sp)]; names(beta.v) <- gsub('beta_', '', x=names(beta.v))
        
    #check growth of nc & meta.rxns combined
    for (thresh in c(0.9, 0.5, 10**-(6:9), 0)){
        rxns1 <- union(names(beta.v)[beta.v>thresh], names(pen.v)[pen.v==0])
        
        cat('Testing if new rxn set gives FBA growth for Sv>=0 w/ threshold', thresh, ': ')
        fba.svg0 <- FBA(sp[,rxns1], fba.obj.rxns=goals, control=list(trace=0), sense='G')
        
        cat('Testing if new rxn set gives FBA growth for Sv=0 w/ threshold', thresh, ': ')
        fba <- FBA(sp[,rxns1], fba.obj.rxns=goals, control=list(trace=0), sense='E')
        cat('\n')
        
        if ((se[1]=='E' & fba$obj>10**-6)|(se[1]=='G' & fba.svg0$obj>10**-6)){
            need.bad.rxns <- intersect(rxns1, names(pen.v)[pen.v>0])
            print(paste("Need bad rxns = c(", paste(need.bad.rxns, collapse="\', \'"), ")", sep="\'"))
            ret <- matrix(cbind(beta.v[need.bad.rxns], pen.v[need.bad.rxns]), ncol=2, dimnames=list(need.bad.rxns, c('beta', 'pen')))
            (ret.o <- ret[order(-ret[,1]),])
            break
        }
    }#end for
    return(ret.o)
}

min.pen.milp <- function(sp, pen.v, model=names(pen.v)[pen.v==0], rxns0=NULL, goals='biomass', f.ub=rep(10**3, ncol(sp)), lb.goals=1, n.add=0, se=rep('E', nrow(sp)), check.input=TRUE, ctrl=list(tilim=60)){
    if (check.input){
        stopifnot(names(pen.v)==colnames(sp), goals %in% colnames(sp), !is.na(pen.v), all(!is.na(sp)), !is.na(f.ub), is.finite(f.ub), model %in% colnames(sp), pen.v[model]==0, is.null(rxns0)|all(rxns0 %in% colnames(sp)))
        cat('Done input validity checks.\n')
    }
    
    #A=[s 0; I -f.ub; 0 gdls.vector]
    #decision vars are [v beta]
    obj <- c(numeric(ncol(sp)), pen.v)
    gdls.v <- numeric(ncol(sp)); names(gdls.v) <- colnames(sp)
    gdls.v[setdiff(names(pen.v)[pen.v>0], rxns0)] <- 1
    
    amat <- rBind(cBind(sp, Matrix(0, nrow=nrow(sp), ncol=ncol(sp))), 
    cBind(Diagonal(n=ncol(sp)), -1*Diagonal(x=f.ub)),
    cBind(Matrix(0, nrow=1, ncol=ncol(sp)), gdls.v))
    
    rhs <- c(numeric(nrow(sp)+ncol(sp)), n.add)
    ub <- rep(c(Inf, 1), each=ncol(sp))
    lb <- numeric(2*ncol(sp))
    #1st sense is for Sv (>)= 0
    sense <- c(se, rep('L', times=ncol(sp)+1))
    #sense <- rep(c(se, 'L'), times=c(nrow(sp), ncol(sp)+1))
    #names
    rownames(amat) <- names(sense) <- names(rhs) <- c(rownames(sp), colnames(sp), 'gdls')
    colnames(amat) <- names(lb) <- names(ub) <- names(obj) <- c(colnames(sp), paste('beta_', colnames(sp), sep=''))
    
    lb[goals] <- lb.goals
    lb[paste('beta_', model, sep='')] <- 1
    vt <- rep(c('C', 'B'), each=ncol(sp))
    
    ##run milp
    mp.milp <- Rcplex(cvec=obj, Amat=amat, bvec=rhs, ub=ub, lb=lb, sense=sense, vtype=vt, control=ctrl)
    cat('CPLEX status of MILP:', mp.milp$stat, '\n')
    x <- mp.milp$xopt; names(x) <- names(obj)
    v <- x[1:ncol(sp)]
    beta.v <- x[ncol(sp) + 1:ncol(sp)]; names(beta.v) <- gsub('beta_', '', x=names(beta.v))
        
    #check growth of nc & meta.rxns combined
    for (thresh in c(0.9, 10**-(6:9), 0)){
        rxns1 <- union(names(beta.v)[beta.v>thresh], names(pen.v)[pen.v==0])
        
        cat('Testing if new rxn set gives FBA growth for Sv>=0 w/ threshold', thresh, ': ')
        fba.svg0 <- FBA(sp[,rxns1], control=list(trace=0), sense='G')
        
        cat('Testing if new rxn set gives FBA growth for Sv=0 w/ threshold', thresh, ': ')
        fba <- FBA(sp[,rxns1], control=list(trace=0), sense='E')
        cat('\n')
        
        if (fba$obj>10**-6){
            need.bad.rxns <- intersect(rxns1, names(pen.v)[pen.v>0])
            print(paste("Need bad rxns = c(", paste(need.bad.rxns, collapse="\', \'"), ")", sep="\'"))
            ret <- matrix(cbind(beta.v[need.bad.rxns], pen.v[need.bad.rxns]), ncol=2, dimnames=list(need.bad.rxns, c('beta', 'pen')))
            (ret.o <- ret[order(-ret[,1], rownames(ret)),])
            break
        }
    }#end for
    return(ret.o)
}
