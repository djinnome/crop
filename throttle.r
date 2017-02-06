##jmd
##5.25.11
##throttle.r

##used in farm_vogel.r
ind <- 1; acc <- 0
length(thr <- which(al.beta<0.5 & al.beta>eps))
fba.ub2 <- fba.ub
fba.ub2[thr] <- fba.ub2[thr]*al.beta[thr]
al.a[vlim.rows, beta.cols] <- -1*Diagonal(x=fba.ub2)
while (acc<0.63 & length(thr)>0 & ind<20){ 
    ##lp
    alarm.lp <- Rcplex(cvec=al.obj, Amat=al.a, bvec=al.rhs, ub=al.ub, lb=al.lb, sense=al.sense, objsense='max')
    stopifnot(alarm.lp$stat==1)
    
    ##extract soln
    x <- alarm.lp$xopt; names(x) <- names(al.obj)
    v <- x[1:nrxns]
    al.beta <- x[beta.cols]; names(al.beta) <- sub('beta_', '', names(al.beta))
    (sum(al.beta<1-eps & al.beta>eps))
    rxns1 <- names(al.beta)[al.beta>eps]
    
    ##test
    ck <- check.ko(s.sp[,rxns1], rxns[rxns1,])
    tab <- table(ck$obs, ck$pred>eps)
    (acc <- tab[1]/sum(tab))
    
    ##throttle bounds
    length(thr <- which(al.beta<1-eps & al.beta>eps))
    fba.ub2[thr] <- fba.ub2[thr]*al.beta[thr]
    al.a[vlim.rows, beta.cols] <- -1*Diagonal(x=fba.ub2)
    
    ind <- ind+1
    flush(stdout())
}

##from farm_seed.r
#get new fba.ub=beta*fba.ub for beta<=0.5 & dual.ub=dual.ub*(1-beta) for beta>0.5
if (!all(c('fva.ub2', 'dual.ub2') %in% ls())){ fva.ub2 <- fva.ub; dual.ub2 <- dual.ub } 
#be conservative with fva.ub to avoid opt_infeas status
fva.ub2[al.beta<=0.5] <- (fva.ub2*(al.beta+10^-3))[al.beta<=0.5]
#can be less conservative with dual, since its objective is not associated with hard constraints
dual.ub2[al.beta>0.5] <- (dual.ub2*(1-al.beta))[al.beta>0.5]
#assign primal ub in al.a
al.a[cbind(nmets + 1:nrxns, beta.cols)] <- -fva.ub2
#assign dual ub in al.a
for (death.cond.i in 1:n.death.conds){
    dci.ko.rows <- nmets+nrxns+2*nrxns*(death.cond.i-1) + nrxns + 1:nrxns
    al.a[cbind(dci.ko.rows, beta.cols)] <- dual.ub2
    al.rhs[dci.ko.rows] <- dual.ub2
}
#run lp
