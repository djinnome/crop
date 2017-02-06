##jmd
##2.8.11
##check_ko_ec.r

check.ko.ec <- function(s.test, rxns.test, rxn.names=rownames(rxns.test), ko.in, sense='G', ub=rep(1000,ncol(s.test)), gapless=FALSE){
    ko <- cbind(ko.in, obs=0, pred=NA)
    
    ##set-up LP
    nmets <- nrow(s.test)
    nrxns <- ncol(s.test)
    obj <- numeric(nrxns)
    obj[rxn.names=='biomass'] <- 1
    names(ub) <- rxn.names
    c.names <- c('CIT-IN-RXN-L2R', 'SUCROSE-IN-RXN-L2R')
    if (all(c.names %in% names(ub))) ub[c.names] <- 50
    
    ##gapless
    if (gapless){ s.test[s.test!=0] <- s.test[s.test!=0]-10^-6 }

    ##loop thru ko
    for (i in 1:nrow(ko)){
        ec.v <- unlist(strsplit(ko$EC[i], split=','))
        rxns.ind.ko <- unique(unlist(apply(as.matrix(ec.v), MARGIN=1, FUN=function(y){ grep(x=rxns.test$EC, pattern=y) })))
        if (length(rxns.ind.ko)>=1){
            ub.ko <- ub
            #KO these rxns by setting ub to 0 (all rxns >= 0 since 2 rxns per reversible one)
            ub.ko[rxns.ind.ko] <- 0
            #sv=0
            cp.tmp <- Rcplex(cvec=obj, Amat=s.test, bvec=rep(0, nrow(s.test)), lb=rep(0, nrxns), ub=ub.ko, objsense='max', sense=sense, control=list(trace=0))
            #check if growth
            if (cp.tmp$status==1){ ko$pred[i] <- cp.tmp$obj } else { print(paste('row', i, 'has cplex status', cp.tmp$status)) }
        }
    }
    return(ko)
}
