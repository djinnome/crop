##jmd
##4.3.11
##alarm_death.r
##apply alarm w/ cutting planes & fva ub's

##PROBLEM: dual is for sv=0 but test is on sv>=0

setwd('/msc/neurospora/FBA/farm')
source('check_ko.r')
source('check_fba.r')
source('cut_mets.r')
source('cut_revs.r')
source('fba.r')
source('farm_header.r')
setwd('/msc/neurospora/FBA/farm_data')

##filter rxns  that can't have flux on vogel
fva.ub0 <- read.csv('s_bal_fva_ub.csv')
fva.ub1 <- fva.ub0[,2]; names(fva.ub1) <- fva.ub0[,1]
rxns.al <- rxns[names(fva.ub1)[fva.ub1>0],]
s.al <- s.sp[,names(fva.ub1)[fva.ub1>0]]
s.al <- s.al[rowSums(abs(s.al))>0,]
#update fva.ub
fva.ub <- fva.ub1[colnames(s.al)]
#update params
nmets <- nrow(s.al); nrxns <- ncol(s.al)

##get experimental KO's
rad <- read.delim('/msc/neurospora/FBA/Neurospora/eCompendium/phenotype-experiment.txt', as.is=TRUE)
rad <- rad[-grep('\\.-', rad$EC),]
rad.rxns.lst <- apply(rad, 1, FUN=function(x){ rxns.al$rxn[rxns.al$EC %in% unlist(strsplit(x['EC'], split=','))] })
rad.rxns.lst <- rad.rxns.lst[-which(sapply(rad.rxns.lst, FUN=function(x) length(x)==0))]
rad.rxns.lst <- rad.rxns.lst[!duplicated(rad.rxns.lst)]

##set list of death conds
#can't yet add sucrose, i think
nut.rxns.lst <- list('OXYGEN-MOLECULE-TRANS-RXN-L2R', 'THIOREDOXIN-TRANS-RXN-L2R', 'BIOTIN-TRANS-RXN-L2R', 
c('AMMONIUM-TRANS-RXN-L2R', 'NITRATE-TRANS-RXN-L2R'))
#set list of rxns
ko.rxns.lst <- c(nut.rxns.lst, rad.rxns.lst)
(n.death.conds <- length(ko.rxns.lst)) #nd

##get cutting planes
cut.met.mat <- cut.mets(s.al)
cut.rev.mat <- cut.revs(s.al, rxns.al)

##instantiate
#decision vars are [v mu1 lambda1 mu2 lambda2 ... beta]
al.a <- Matrix(0, nrow=nmets+nrxns+2*nrxns*n.death.conds+nrow(cut.met.mat)+nrow(cut.rev.mat), 
ncol=nrxns+(nmets+nrxns)*n.death.conds+nrxns)
al.lb <- al.obj <- rep(0, ncol(al.a))
al.ub <- rep(Inf, ncol(al.a))
al.rhs <- rep(0, nrow(al.a))
#set sense as NA and fill in below
al.sense <- character(nrow(al.a))
#index beta
beta.cols <- nrxns+n.death.conds*(nmets+nrxns) + 1:nrxns
#set carbon source rnxs
c.rxns <- c('SUCROSE-TRANS-RXN-L2R','CIT-TRANS-RXN-L2R')

##names
#rownames
cut.rownames <- c(paste('cut_met', 1:nrow(cut.met.mat), sep=''), paste('cut_rev', 1:nrow(cut.rev.mat), sep=''))
if (n.death.conds>0){ 
    dual.rownames <- paste('dual_', rep(c('eq', 'ko'), each=nrxns), rep(1:n.death.conds, each=2*nrxns), '_', colnames(s.al), sep='')
    rownames(al.a) <- c(paste('fba', c(rownames(s.al), colnames(s.al)), sep='_'), dual.rownames, cut.rownames)
} else { 
    rownames(al.a) <- c(paste('fba', c(rownames(s.al), colnames(s.al)), sep='_'), cut.rownames)
}
#colnames
if (n.death.conds>0){ 
    dual.colnames <- paste(rep(c('mu','l'), times=c(nmets, nrxns)), rep(1:n.death.conds, each=nmets+nrxns), '_', c(rownames(s.al), colnames(s.al)), sep='')
    colnames(al.a) <- c(colnames(s.al), dual.colnames, paste('beta', colnames(s.al), sep='_'))
} else { 
    colnames(al.a) <- c(colnames(s.al), paste('beta', colnames(s.al), sep='_')) 
}
#vector names
names(al.lb) <- names(al.ub) <- names(al.obj) <- colnames(al.a)
names(al.rhs) <- names(al.sense) <- rownames(al.a)

##growth on vogel
#sv>=0
#by default, sv rhs already 0
al.a[1:nmets, 1:nrxns] <- s.al
al.sense[1:nmets] <- 'G'
#Iv-fva.ub*beta<=0
al.a[nmets + 1:nrxns, 1:nrxns] <- Diagonal(n=ncol(s.al))
al.a[nmets + 1:nrxns, beta.cols] <- -1*Diagonal(x=fva.ub)
#accidentally tried Iv=1000*beta, w/ good results
al.sense[nmets + 1:nrxns] <- 'L'
#bounds on fluxes
#if this ub=100, growth ~= 3
al.ub[c.rxns] <- 10
#bound other nut fluzes st 1000 isn't too low of an ub for other fluxes
al.ub[setdiff(vogel$rxn, c.rxns)] <- 100
al.lb['biomass'] <- 0.1

##run fba, for use in dual
fba.obj <- rep(0, nrxns)
names(fba.obj) <- colnames(s.al)
fba.obj['biomass'] <- 1
#lp
fba.lp <- Rcplex(cvec=fba.obj, Amat=s.al, bvec=rep(0, nrow(s.al)), ub=fva.ub, lb=0, sense='G', objsense='max')
fba.lp$stat; fba.obj.val <- fba.lp$obj
#dual ub: check w/ paschalidis
dual.ub <- fba.obj.val/fva.ub

##beta bounds & coefficients
al.obj[beta.cols] <- 0.5-rxns.al$prob
#by default beta lb already 0
al.ub[beta.cols] <- 1
#set betas of desired fluxes to 1
keep.rxns <- c('biomass', vogel$rxn, unlist(ko.rxns.lst))
al.lb[paste('beta', unique(keep.rxns), sep='_')] <- 1

##duals
#much faster to define in loop instead of defining in fba.dual & passing
#loop
for (death.cond.i in 1:n.death.conds){
    ko.rxns <- ko.rxns.lst[[death.cond.i]]
    
    ##get indices
    #rows for eq: s'mu-lambda=-c & ko: lambda+dual.ub*beta<=dual.ub
    dci.rows <- nmets+nrxns+2*nrxns*(death.cond.i-1) + 1:(2*nrxns)
    dci.eq.rows <- nmets+nrxns+2*nrxns*(death.cond.i-1) + 1:nrxns
    dci.ko.rows <- nmets+nrxns+2*nrxns*(death.cond.i-1) + nrxns + 1:nrxns
    #cols for dual vars mu & lambda
    dci.cols <- nrxns+(death.cond.i-1)*(nmets+nrxns) + 1:(nmets+nrxns)
    dci.mu.cols <- nrxns+(death.cond.i-1)*(nmets+nrxns) + 1:nmets
    dci.l.cols <- nrxns+(death.cond.i-1)*(nmets+nrxns) + nmets + 1:nrxns
    
    ##assign dual submatrices
    #eq & ko rows same for all dci.conds but need to break into mult assignments or else really slow
    al.a[dci.eq.rows, dci.mu.cols] <- t(s.al)
    al.a[dci.eq.rows, dci.l.cols] <- -1*Diagonal(ncol(s.al))
    al.a[dci.ko.rows, dci.l.cols] <- Diagonal(nrxns)
    
    ##assign betas
    dci.beta.coeff <- Diagonal(x=dual.ub)
    #don't want beta=1 -> lambda<0 for ko.rxns or nutrient rxns
    rxns.no.beta <- c(vogel$rxn, ko.rxns)
    diag(dci.beta.coeff)[which(colnames(s.al) %in% rxns.no.beta)] <- 0
    al.a[dci.ko.rows, beta.cols] <- dci.beta.coeff
    
    ##bounds on duals
    #ub=Inf by default, above
    al.lb[dci.l.cols] <- -Inf
    #mu>=0 for sv>=0
    al.lb[dci.mu.cols] <- 0
            
    ##bounds & penalties on nut rxns
    al.lb[paste('l', death.cond.i, '_', vogel$rxn, sep='')] <- 0
    al.obj[paste('l', death.cond.i, '_', vogel$rxn, sep='')] <- 100
    al.obj[paste('l', death.cond.i, '_', c.rxns, sep='')] <- 1000
    
    ##rhs
    al.sense[dci.rows] <- rep(c('E', 'L'), each=nrxns)
    al.rhs[dci.rows] <- c(-fba.obj, dual.ub)
}

##cutting plane matrices
al.a[nmets+nrxns+2*nrxns*n.death.conds + 1:nrow(cut.met.mat), beta.cols] <- cut.met.mat
al.a[nmets+nrxns+2*nrxns*n.death.conds+nrow(cut.met.mat)+ 1:nrow(cut.rev.mat), beta.cols] <- cut.rev.mat
al.sense[nmets+nrxns+2*nrxns*n.death.conds + 1:(nrow(cut.met.mat)+nrow(cut.rev.mat))] <- 'L'

##alarm lp
#min weighted penalty of beta's
alarm.lp <- Rcplex(cvec=al.obj, Amat=al.a, bvec=al.rhs, ub=al.ub, lb=al.lb, sense=al.sense)
alarm.lp$stat

##extract soln
#x
x <- alarm.lp  $xopt; names(x) <- names(al.obj)
#beta
al.beta <- x[(length(x)-nrxns+1):length(x)]
eps <- 10^-6
mean(al.beta>1-eps|al.beta<eps)
sum(al.beta>eps)
rxns1 <- sub('beta_', '', names(al.beta)[al.beta>eps])

##vogel growth
vv <- fba(s.al[,rxns1], sense='G')

##throttle
#paste these bounds into al.a (except don't run dual.ub=fba.growth.obj/fva.ub) then run alarm.lp
#primal
fva.ub.ind <- which(al.beta>=0.5 & al.beta<1-eps)
fva.ub[fva.ub.ind] <- al.beta[fva.ub.ind]*fva.ub[fva.ub.ind]
#dual
dual.ub.ind <- which(al.beta<0.5 & al.beta>eps)
dual.ub[dual.ub.ind] <- al.beta[dual.ub.ind]*dual.ub[dual.ub.ind]

##death.conds
#beta=0.002 becomes beta=1 (n crassa can't grow on vogel's if set these to 0) so growth jumps dramatically from lambda*beta
#vogel lambda's
x[paste('l', 1:n.death.conds, '_', vogel$rxn, sep='')]
#lambda*fva.ub*beta
lam <- x[paste('l', death.cond.i, '_', colnames(s.al), sep='')]
prd <- lam*fva.ub*al.beta; names(prd) <- colnames(s.al)
#only want product from v<=fva.ub*beta, which have lambda>=0
prd <- prd[prd>0]; prd[order(-prd)]
#check fba's
for (death.cond.i in 1:n.death.conds){
    ko.rxns <- ko.rxns.lst[[death.cond.i]]
    if (all(ko.rxns %in% rxns1)){
        fba.tmp <- fba(s.al[,rxns1], ko.rxns=ko.rxns, sense='G', fba.ub=fva.ub[rxns1])
        cat('KOs', ko.rxns, 'status', fba.tmp$status, 'w/ fba.obj', fba.tmp$obj, '\n \n')
    }
}

##test
#check fba
(pred <- check.fba(a=s.al[,rxns1], ko.rxns=ko.rxns.lst, sense='G'))
mean(pred<0.1)
#check ko
ck <- check.ko(s.al[,rxns1], rxns.al[rxns1,])
table(ck$obs, as.numeric(ck$pred>eps))

##write
rr.out <- data.frame(frame=rxns.al[rxns1, 3], flux=x[rxns1], rxns.al[rxns1, -3])
write.table(rr.out, 'rxns_al_death.txt', quote=FALSE, row.names=FALSE, sep='\t')

###validate##############################################################################################################
##dual in loop
col.inds <- c(dci.cols, beta.cols)
dual.lp <- Rcplex(cvec=c(al.obj[dci.cols], rep(0, nrxns)), Amat=al.a[dci.rows, col.inds], 
bvec=dual.rhs, sense=dual.sense, lb=c(al.lb[dci.cols], rep(1, length(beta.cols))), ub=al.ub[col.inds])
dual.lp$stat
dual.lp$obj
