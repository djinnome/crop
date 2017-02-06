##jmd
##3.8.11
##alarm_script.r
##this version for universal matrix w/ growth in 1 cond and no imports for non-media

#call: bsub -q week -P nc10 -o foo.log /home/radon00/jdreyf/farm/alarm_script.r

setwd('/home/radon00/jdreyf/farm')

source('farm_header.r')

##get experimental KO's
rad <- read.delim('/msc/neurospora/FBA/Neurospora/eCompendium/phenotype-experiment.txt', as.is=TRUE)
rad <- rad[-grep('\\.-', rad$EC),]
rad.rxns.lst <- apply(rad, 1, FUN=function(x){ rxns$rxn[rxns$EC %in% unlist(strsplit(x['EC'], split=','))] })
rad.rxns.lst <- rad.rxns.lst[-which(sapply(rad.rxns.lst, FUN=function(x) length(x)==0), arr.ind=TRUE)]
rad.rxns.lst <- rad.rxns.lst[!duplicated(rad.rxns.lst)]

##set list of death conds
#can't yet add sucrose, i think
ko.rxns.lst <- c(list('OXYGEN-MOLECULE-TRANS-RXN-L2R', 'THIOREDOXIN-TRANS-RXN-L2R', 'BIOTIN-TRANS-RXN-L2R', 
c('AMMONIUM-TRANS-RXN-L2R', 'NITRATE-TRANS-RXN-L2R')), rad.rxns.lst)
(n.death.conds <- length(ko.rxns.lst)) #nd
# stopifnot(length(ko.rxns.lst)==n.death.conds)

##decision vars are [v mu1 lambda1 mu2 lambda2 beta] where 2 dual condns are w/o thioredoxin & o2
#obj: min (1-p)beta+1000*(l1_cit+l1_suc+l2_cit+l2_suc)
#constraints: cv>=2; sv=0; v-M*beta<=0
#dual.amat1=rBind(eq,ko.mat1) (sense) dual.rhs
#dual.amat2=rBind(eq,ko.mat2) (sense) dual.rhs

##1 life cond & nd>2 death conds
#obj: min [0_(n+(m+n)*nd) (1-p)_n] #except l(i)_(cit|suc)=10^3
#fba
#[1 row: c_(m) 0_(nd*(m+n)+n) >= 2;
#m rows: s_(n) 0_(nd*(m+n)+n) = 0;
#death cond i
#n rows: 0_(n+(i-1)*(m+n)) t(s)_(m) -1*I_(n) 0_((m+n)*(nd-i)+n) = -c; #eq
#m rows: 0_(n+(n+m)*(i-1)+m) I_(n) 0_((m+n)*(nd-i)) 1000*I_(n) <= 1000; #ko.mat: except rows in {carbon src, a priori ko}

##instantiate
# m <- nmets; n <- nrxns
al.a <- Matrix(0, nrow=1+nmets+nrxns+2*nrxns*n.death.conds, ncol=nrxns+(nrxns+nmets)*n.death.conds+nrxns)
al.lb <- al.obj <- rep(0, ncol(al.a))
al.ub <- rep(Inf, ncol(al.a))
al.rhs <- rep(0, nrow(al.a))
al.sense <- rep('E', nrow(al.a))
#index beta
beta.cols <- nrxns+n.death.conds*(nmets+nrxns) + 1:nrxns

##names
#rownames
if (n.death.conds>0){ 
    dual.rownames <- paste('dual_', rep(c('eq', 'ko'), each=nrxns), rep(1:n.death.conds, each=2*nrxns), '_', colnames(s.sp), sep='')
    rownames(al.a) <- c('primal.obj', paste('fba', c(rownames(s.sp), colnames(s.sp)), sep='_'), dual.rownames)
} else { 
    rownames(al.a) <- c('primal.obj', paste('fba', c(rownames(s.sp), colnames(s.sp)), sep='_'))
}
#colnames
dual.colnames <- paste(rep(c('mu','l'), times=c(nmets, nrxns)), rep(1:n.death.conds, each=nmets+nrxns), '_', c(rownames(s.sp), colnames(s.sp)), sep='')
if (n.death.conds>0){ 
    colnames(al.a) <- c(colnames(s.sp), dual.colnames, paste('beta', colnames(s.sp), sep='_'))
} else { 
    colnames(al.a) <- c(colnames(s.sp), paste('beta', colnames(s.sp), sep='_')) 
}
#vector names
names(al.lb) <- names(al.ub) <- names(al.obj) <- colnames(al.a)
names(al.rhs) <- names(al.sense) <- rownames(al.a)

##run fba
c.rxns <- c('SUCROSE-TRANS-RXN-L2R','CIT-TRANS-RXN-L2R')
fba.res <- fba(a=s.sp, c.rxns=c.rxns, fba.run=FALSE)

#fba
al.a[1,1:nrxns] <- fba.res$obj.v
al.a[2:(1+nmets), 1:nrxns] <- s.sp
al.a[(2+nmets):(1+nmets+nrxns), 1:nrxns] <- Diagonal(ncol(s.sp))
al.a[(2+nmets):(1+nmets+nrxns), (ncol(al.a)-nrxns+1):ncol(al.a)] <- -1000*Diagonal(ncol(s.sp))
#bounds
al.rhs[1] <- 2.5
al.sense[1] <- 'G'
al.ub[1:nrxns] <- fba.res$ub

##duals
#much faster to define in loop instead of defining in fba.dual & passing
#set vars which are invariant to cond'n
dual.sense <- rep(c('E', 'L'), each=nrxns)
dual.rhs <- c(-fba.res$obj.v, rep(10^6, nrxns))
#loop
for (death.cond.i in 1:n.death.conds){
    ko.rxns <- ko.rxns.lst[[death.cond.i]]
    
    ##get indices
    #rows for eq: s'mu-l=-c & ko: l+1000*beta<=1000
    dci.rows <- 1+nmets+nrxns+2*nrxns*(death.cond.i-1) + 1:(2*nrxns)
    dci.eq.rows <- 1+nmets+nrxns+2*nrxns*(death.cond.i-1) + 1:nrxns
    dci.ko.rows <- 1+nmets+nrxns+2*nrxns*(death.cond.i-1) + nrxns + 1:nrxns
    #cols for dual vars mu & lambda
    dci.cols <- nrxns+(death.cond.i-1)*(nmets+nrxns) + 1:(nmets+nrxns)
    dci.mu.cols <- nrxns+(death.cond.i-1)*(nmets+nrxns) + 1:nmets
    dci.l.cols <- nrxns+(death.cond.i-1)*(nmets+nrxns) + nmets + 1:nrxns
    
    ##assign dual submatrices
    #eq & ko rows same for all dci.cols but need to break into mult assignments or else really slow!
    al.a[dci.eq.rows, dci.mu.cols] <- t(s.sp)
    al.a[dci.eq.rows, dci.l.cols] <- -1*Diagonal(ncol(s.sp))
    al.a[dci.ko.rows, dci.l.cols] <- Diagonal(nrxns)
    
    ##assign betas
    dci.beta.coeff <- 1000*Diagonal(nrxns)
    diag(dci.beta.coeff)[which(colnames(s.sp) %in% c(c.rxns, ko.rxns))] <- 0
    al.a[dci.ko.rows, beta.cols] <- dci.beta.coeff
    
    ##bounds on duals
    al.lb[dci.cols] <- -Inf
    al.ub[dci.cols] <- Inf
    #c.rxns
    al.lb[paste('l', death.cond.i, '_', c.rxns, sep='')] <- 0
    al.obj[paste('l', death.cond.i, '_', c.rxns, sep='')] <- 1000
    #rhs, from just before loop
    al.sense[dci.rows] <- dual.sense
    al.rhs[dci.rows] <- dual.rhs
}

##beta's
al.obj[beta.cols] <- 0.99-rxns$prob
al.lb[beta.cols] <- 0
al.ub[beta.cols] <- 1
#set desired fluxes to 1
# al.lb[paste('beta', vogel$rxn, sep='_')] <- 1

##alarm lp
alarm.lp <- Rcplex(cvec=al.obj, Amat=al.a, bvec=al.rhs, ub=al.ub, lb=al.lb, sense=al.sense)

#write
rr.out <- data.frame(frame=rxns[rxns1, 3], flux=x[rxns1], rxns[rxns1, -3])
write.table(rr.out, 'rxns_beta1.txt', quote=FALSE, row.names=FALSE, sep='\t')