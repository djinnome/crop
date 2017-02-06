##jmd
##6.5.11
##farm_seed.r
##apply farm to seed db, uses sv=0

library('Rcplex')
library('Matrix')

setwd('/msc/neurospora/FBA/farm')
source('check_ko.r')
source('check_fba.r')
source('cut_mets.r')
source('cut_revs.r')
source('fba.r')
setwd('../seed')
options(stringsAsFactors=FALSE)

##read
bm <- read.delim('biomass_seed.txt')
rxns <- read.csv('rxns_seed_annot.csv', row.names=1)
vogel <- read.delim('vogel_seed.txt')
s <- read.delim('s_seed.txt', as.is=FALSE)
s.sp <- sparseMatrix(i=as.numeric(s$compound), j=as.numeric(s$rxn), x=as.real(s$coeff))
dimnames(s.sp) <- list(levels(s$compound), levels(s$rxn))
nmets <- nrow(s.sp); nrxns <- ncol(s.sp)
#fva ub
fva.ub <- rep(1000, ncol(s.sp))
names(fva.ub) <- colnames(s.sp)
#rename, for below
s.al <- s.sp
rxns.al <- rxns

##get experimental KO's
rad <- read.delim('/msc/neurospora/FBA/Neurospora/eCompendium/phenotype-experiment.txt', as.is=TRUE)
rad <- rad[!duplicated(rad$EC),]
#get rid of 3-level ec
rad <- rad[-grep('\\.-', rad$EC),]
#make rxn list
rad.rxns.lst <- list()
for (i in 1:nrow(rad)){
    ec.v <- unlist(strsplit(rad$EC[i], split=','))
    rxns.ind.ko <- unique(unlist(apply(as.matrix(ec.v), MARGIN=1, FUN=function(y){ rownames(rxns.al)[grep(x=rxns.al$EC, pattern=y)] })))
    if (length(rxns.ind.ko)>=1){ rad.rxns.lst[[i]] <- as.character(rxns.ind.ko) }
}
#rad.rxns.lst <- rad.rxns.lst[-which(sapply(rad.rxns.lst, FUN=function(x) length(x)==0))]
rad.rxns.lst <- rad.rxns.lst[!duplicated(rad.rxns.lst)]

##set list of death conds
#can't yet add sucrose, i think
nut.rxns.lst <- list('OXYGEN-MOLECULE-TRANS-RXN-L2R')
#rad.rxns.lst <-  rad.rxns.lst[1]
rad.rxns.lst <- NULL
#set list of rxns
ko.rxns.lst <- c(nut.rxns.lst, rad.rxns.lst)
(n.death.conds <- length(ko.rxns.lst))

##get cutting planes
cut.met.mat <- cut.mets(s.sp, sense='E')
#cut.met.mat <- Matrix(0, nrow=1, ncol=nrxns)
#cut.rev.mat <- cut.revs(s.sp, rxns)
cut.rev.mat <- Matrix(0, nrow=1, ncol=nrxns)

##instantiate
#conds: vogel & gapless growth, and mult death
#growth cond has nmets+nrxns rows, for sv=0 & v<=M*beta
#decision vars are [v1 mu1 lambda1 mu2 lambda2 ... v2 beta]
#constraints are [vogel.growth death.conds gapless cut.planes]
al.a <- Matrix(0, nrow=nmets+nrxns+2*nrxns*n.death.conds+nmets+nrxns+nrow(cut.met.mat)+nrow(cut.rev.mat), 
ncol=nrxns+(nmets+nrxns)*n.death.conds+nrxns+nrxns)
al.lb <- al.obj <- rep(0, ncol(al.a))
al.ub <- rep(Inf, ncol(al.a))
al.rhs <- rep(0, nrow(al.a))
#set sense as NA and fill in below
al.sense <- character(nrow(al.a))
#index beta
beta.cols <- 2*nrxns+n.death.conds*(nmets+nrxns) + 1:nrxns
#set carbon source rnxs
c.rxns <- c('SUCROSE-TRANS-RXN-L2R','CIT-TRANS-RXN-L2R')

##names
#rownames
gapless.rownames <- paste('gapless', c(rownames(s.al), colnames(s.al)), sep='_')
cut.rownames <- c(paste('cut_met', 1:nrow(cut.met.mat), sep=''), paste('cut_rev', 1:nrow(cut.rev.mat), sep=''))
if (n.death.conds>0){ 
    dual.rownames <- paste('dual_', rep(c('eq', 'ko'), each=nrxns), rep(1:n.death.conds, each=2*nrxns), '_', colnames(s.al), sep='')
    rownames(al.a) <- c(paste('fba', c(rownames(s.al), colnames(s.al)), sep='_'), dual.rownames, gapless.rownames, cut.rownames)
} else { 
    rownames(al.a) <- c(paste('fba', c(rownames(s.al), colnames(s.al)), sep='_'), gapless.rownames, cut.rownames)
}
#colnames
if (n.death.conds>0){ 
    dual.colnames <- paste(rep(c('mu','l'), times=c(nmets, nrxns)), rep(1:n.death.conds, each=nmets+nrxns), '_', c(rownames(s.al), colnames(s.al)), sep='')
    colnames(al.a) <- c(colnames(s.al), dual.colnames, paste('gapless', colnames(s.al), sep='_'), paste('beta', colnames(s.al), sep='_'))
} else { 
    colnames(al.a) <- c(colnames(s.al), paste('beta', colnames(s.al), sep='_')) 
}
#vector names
names(al.lb) <- names(al.ub) <- names(al.obj) <- colnames(al.a)
names(al.rhs) <- names(al.sense) <- rownames(al.a)

##growth on vogel
#sv>=0
#by default, sv rhs already 0
s.gapless <- s.al
s.gapless[s.gapless!=0] <- s.gapless[s.gapless!=0]-10^(-3)
al.a[1:nmets, 1:nrxns] <- s.gapless
#make 'E' for sv=0
al.sense[1:nmets] <- 'E'
#Iv-fva.ub*beta<=0
al.a[nmets + 1:nrxns, 1:nrxns] <- Diagonal(n=ncol(s.al))
al.a[nmets + 1:nrxns, beta.cols] <- -1*Diagonal(x=fva.ub)
#accidentally tried Iv=1000*beta, w/ good results
al.sense[nmets + 1:nrxns] <- 'L'
#bounds on fluxes
#if this ub=100, growth ~= 3
al.ub[c.rxns] <- 1000
#bound other nut fluzes st 1000 isn't too low of an ub for other fluxes
al.ub[setdiff(vogel$rxn, c.rxns)] <- 1000
al.lb['biomass'] <- 0.05

##run fba, for use in dual
fba.obj <- rep(0, nrxns)
names(fba.obj) <- colnames(s.al)
fba.obj['biomass'] <- 1
#lp
fba.lp <- Rcplex(cvec=fba.obj, Amat=s.al, bvec=rep(0, nrow(s.al)), ub=fva.ub, lb=0, sense='G', objsense='max')
fba.lp$stat
(fba.obj.val <- fba.lp$obj)
#dual ub: check w/ paschalidis
dual.ub <- fba.obj.val/fva.ub

##beta bounds & coefficients
rxns.al[unlist(ko.rxns.lst), 'prob'] <- 1
al.obj[beta.cols] <- 0.95-rxns.al$prob
#by default beta lb already 0
al.ub[beta.cols] <- 1
#set betas of desired fluxes to 1
#keep.rxns <- c('biomass', vogel$rxn, unlist(ko.rxns.lst))
#al.lb[paste('beta', c.rxns, sep='_')] <- 1

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
    #mu free for sv=0
    al.lb[dci.mu.cols] <- -Inf
            
    ##bounds & penalties on nut rxns
    al.lb[paste('l', death.cond.i, '_', vogel$rxn, sep='')] <- 0
    al.obj[paste('l', death.cond.i, '_', vogel$rxn, sep='')] <- 100
    al.obj[paste('l', death.cond.i, '_', c.rxns, sep='')] <- 1000
    
    ##rhs
    al.sense[dci.rows] <- rep(c('E', 'L'), each=nrxns)
    al.rhs[dci.rows] <- c(-fba.obj, dual.ub)
    
    print(death.cond.i)
}

##gapless
#Sv-eps*|S|*beta>=0
#v-fva.ub*beta<=0
gl.v.cols <- nrxns+(nmets+nrxns)*n.death.conds + 1:nrxns
gl.produce.rows <- nmets+nrxns+2*nrxns*n.death.conds + 1:nmets
gl.ub.rows <- nmets+nrxns+2*nrxns*n.death.conds+nmets + 1:nrxns
#assign
al.a[gl.produce.rows, gl.v.cols] <- sign(s.al)
al.a[gl.produce.rows, beta.cols] <- 10^(-3)*abs(sign(s.al))
al.a[gl.ub.rows, gl.v.cols] <- Diagonal(n=ncol(s.al))
al.a[gl.ub.rows, beta.cols] <- -1*Diagonal(x=fva.ub)
#rhs
#need to make all nuts highly available to capture feeder pathways for nuts not in vogel
#don't currently have rxns bringing these nuts in extracellular compartment, so allow these extracellular mets to deplete
al.rhs[gl.produce.rows][grep('[e]', names(al.rhs[gl.produce.rows]), fixed=TRUE)] <- -1000
#sense
al.sense[gl.produce.rows] <- 'G'
al.sense[gl.ub.rows] <- 'L'

##cutting plane matrices
al.a[nmets+nrxns+2*nrxns*n.death.conds+nmets+nrxns + 1:nrow(cut.met.mat), beta.cols] <- cut.met.mat
al.a[nmets+nrxns+2*nrxns*n.death.conds+nmets+nrxns+nrow(cut.met.mat) + 1:nrow(cut.rev.mat), beta.cols] <- cut.rev.mat
al.sense[nmets+nrxns+2*nrxns*n.death.conds+nmets+nrxns + 1:(nrow(cut.met.mat)+nrow(cut.rev.mat))] <- 'L'

##alarm lp
#min weighted penalty of beta's
alarm.lp <- Rcplex(cvec=al.obj, Amat=al.a, bvec=al.rhs, ub=al.ub, lb=al.lb, sense=al.sense)
alarm.lp$stat

##extract soln
#x
x <- alarm.lp$xopt; names(x) <- names(al.obj)
v <- x[1:nrxns]
#beta
al.beta <- x[(length(x)-nrxns+1):length(x)]
eps <- 10^-6
mean(al.beta>1-eps|al.beta<eps)
sum(al.beta>eps)
rxns1 <- sub('beta_', '', names(al.beta)[al.beta>0|v>0])
#inconsistent flux/beta
v[v>0 & al.beta==0]

##vogel growth
fba.r1.ub=rep(10^3, length(rxns1)); names(fba.r1.ub) <- rxns1
vv <- FBA(s.al[,rxns1], sense='E', fba.ub=fba.r1.ub)
#vv <- FBA(s.al[,rxns1], sense='E', goal.rxns=unique(bm$rxn))
#gapless
vs <- FBA(a=s.gapless[,rxns1], fba.ub=rep(10^3, length(rxns1)), sense='E')
#no vogel
fba.r1.ub[names(fba.r1.ub) %in% vogel$rxn] <- 0
vn <- FBA(a=s.gapless[,rxns1], fba.ub=fba.r1.ub, sense='E')
v <- vn$xopt; names(v) <- rxns1; rxns[intersect(rxns1[v>999], rownames(rxns)[rxns$prob<0]),]

##check ko
#gapless
ck <- check.ko.ec(s.gapless[,rxns1], rxns[rxns1,], ko=rad, sense='E', ub=rep(10^3, length(rxns1)))
table(ck$obs, ck$pred>eps)
#usual
ck <- check.ko.ec(s.sp[,rxns1], rxns[rxns1,], ko=rad)
table(ck$obs, ck$pred>eps)
##write
rr.out <- data.frame(frame=rxns.al[rxns1,'meta'], flux=x[rxns1], beta=al.beta[paste('beta', rxns1, sep='_')], coeff=al.obj[paste('beta', rxns1, sep='_')], row.name=rownames(rxns.al[rxns1,]), rxns.al[rxns1, -3])
#write.table(rr.out, 'rxns_june13.txt', quote=FALSE, row.names=FALSE, sep='\t')

###validate##############################################################################################################
##dual in loop
col.inds <- c(dci.cols, beta.cols)
dual.lp <- Rcplex(cvec=c(al.obj[dci.cols], rep(0, nrxns)), Amat=al.a[dci.rows, col.inds], 
bvec=dual.rhs, sense=dual.sense, lb=c(al.lb[dci.cols], rep(1, length(beta.cols))), ub=al.ub[col.inds])
dual.lp$stat
dual.lp$obj
