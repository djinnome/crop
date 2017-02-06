##jmd
##3.7.11
##fba.r

require('Matrix')
require('Rcplex')

##after farm, some mutants can grow w/ sv>=0 but not sv=0
#a needs to have colnames
#eps is for 'secondary' flux min
#min2 is flag to perform secondary flux min
#ngam dominates fab.lb
FBA <- function(a, fba.obj.rxns='biomass', goals=NULL, ko.rxns=NULL, sense='E', fba.ub=NULL, fba.lb=NULL, min2=FALSE, ngam.name='NGAM', ngam.val=0, 
fbawmc=NULL, eps=0, ub.goals=1, quiet=FALSE, cprim0=TRUE, ...){
    stopifnot(is.null(fba.obj.rxns)|fba.obj.rxns %in% colnames(a), goals %in% c(colnames(a), rownames(a)))
    
    if (is.null(fba.ub)){ fba.ub <- rep(1000, ncol(a)) }
    if (is.null(fba.lb)){ fba.lb <- numeric(ncol(a)) }
    
    if (!all(ko.rxns %in% colnames(a))){ 
        cat('ko rxns:', ko.rxns[!(ko.rxns %in% colnames(a))], 'not in S matrix \n')
        ko.rxns <- ko.rxns[ko.rxns %in% colnames(a)]
        if (length(ko.rxns)==0) ko.rxns <- NULL
    }
    if (!is.null(goals)) fba.obj.rxns  <- NULL
    
    ##goals
    #given goal rxns
    goal.rxns <- goals[goals %in% colnames(a)]
    #construct metabolite goal rxns
    if (any(goals %in% rownames(a))){
        goal.mets <- goals[goals %in% rownames(a)]
        goal.met.mat <- apply(as.matrix(goal.mets), 1, FUN=function(x){ y <- numeric(nrow(a)); names(y) <- rownames(a); y[x] <- -1; return(y) })
        colnames(goal.met.mat) <- paste(goal.mets, 'GOAL-RXN-L2R', sep='-')
        a <- cBind(a, goal.met.mat)
        goal.rxns <- c(goal.rxns, colnames(goal.met.mat))
        
        #cat('ncol a:', ncol(a), 'l ub:', length(fba.ub), '\n')
        #extend bounds to inc. new rxns
        #maybe fba.ub definition in fcn call evaluates below, so that length(fba.ub) = ncol(new a)
        fba.ub <- c(fba.ub, rep(median(fba.ub), length(goal.mets)))
        fba.lb <- c(fba.lb, rep(median(fba.lb), length(goal.mets)))
    }
    #cat('ncol a:', ncol(a), 'l ub:', length(fba.ub), '\n')
    
    fba.obj <- rep(0, ncol(a))
    if (cprim0){ cprim <- numeric(ncol(a)) } else { cprim=NULL }
    names(fba.obj) <- names(fba.lb) <- names(fba.ub) <- colnames(a)
    if (ngam.name %in% colnames(a)){ fba.lb[ngam.name] <- ngam.val }
    
    #fba.ub[intersect(vogel$rxn, names(fba.ub))] <- apply(as.matrix(fba.ub[intersect(vogel$rxn, names(fba.ub))]), 1, FUN=function(x){ min(x, 1000) })
    
    fba.obj[union(fba.obj.rxns, goal.rxns)] <- 1
    #fba.ub[fba.obj.rxns] <- 10^3
    fba.ub[ko.rxns] <- 0
    #goal programming
    fba.ub[goal.rxns] <- ub.goals
    #secondary flux min
    fba.obj[!(names(fba.obj) %in% union(fba.obj.rxns, goal.rxns))] <- -eps
    
    #lb=0 by default
    #fbawmc
    se <- character(nrow(a)); se[1:nrow(a)] <- sense
    if (!is.null(fbawmc)){
        a1 <- rBind(a, rep(1, ncol(a)))
        rownames(a1)[nrow(a1)] <- 'fbawmc'
        se1 <- c(se, 'L')
        
        max.bm <- Rcplex(cvec=fba.obj, Amat=a1, bvec=c(rep(0, nrow(a)), fbawmc), ub=fba.ub, lb=fba.lb, sense=se1, objsense='max', cprim=cprim, ...)
        names(max.bm$xopt) <- colnames(a1)
    } else {
        max.bm <- Rcplex(cvec=fba.obj, Amat=a, bvec=rep(0, nrow(a)), ub=fba.ub, lb=fba.lb, sense=se, objsense='max', cprim=cprim, ...)
        names(max.bm$xopt) <- colnames(a)
    }
    
    if (!(max.bm$stat %in% 1:2)){ max.bm$obj <- NA }   
    ret <- max.bm
    
    if (!quiet){ 
        cat('CPLEX status:', max.bm$status, 'w/ unpenalized objective:', sum(max.bm$xopt[union(fba.obj.rxns, goal.rxns)]), '\n')
        if (length(goals)>=1){ cat('Unmet goals are:', goals[max.bm$xopt[goal.rxns]<10**-6], '\n') }
    }#end if !quiet
    
    ##2nd min
    #add new obj >= previous obj. could make equal but this could have better numerical properties.
    if (!is.na(max.bm$obj) & min2){
        a2 <- rBind(a, fba.obj)
        b2 <- c(rep(0, nrow(a)), max.bm$obj)
        #handle sense whether it's a single character or a vector
        se2 <- character(nrow(a2)); se2[1:nrow(a)] <- sense; se2[nrow(a)+1] <- 'E'
        min.v <- Rcplex(cvec=rep(1,ncol(a)), Amat=a2, bvec=b2, ub=fba.ub, lb=fba.lb, sense=se2, objsense='min', control=list(trace=0))
        names(min.v$xopt) <- colnames(a)
        if (!quiet){ cat('CPLEX status of 2nd min:', min.v$status, 'w/ objective:', min.v$obj, '\n') }
        ret <- min.v
    }
    
    return(ret)
}#end fcn
