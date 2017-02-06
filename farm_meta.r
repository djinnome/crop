##jmd
##6.5.11
##farm_meta.r
##apply farm to meta db, uses sv>=0

library('Rcplex')
library('Matrix')

setwd('/msc/neurospora/FBA/farm')
source('check_ko.r'); source('check_fba.r')
source('cut_mets.r'); source('cut_revs.r')
source('fba.r'); source('fba_dual.r'); source('jdc.r')
source('phenos2rxns.r')
setwd('../farm_data')
options(stringsAsFactors=FALSE)

##read
bm <- read.delim('biomass_s.txt')
vogel <- read.delim('vogel_trans_s.txt')
s <- read.delim('s_bal.txt', as.is=FALSE, na='NIL')
s.sp <- sparseMatrix(i=as.numeric(s$compound), j=as.numeric(s$rxn), x=as.real(s$coeff))
dimnames(s.sp) <- list(levels(s$compound), levels(s$rxn))
nmets <- nrow(s.sp); nrxns <- ncol(s.sp)
#rxns annot
rxns <- read.delim('rxn_bal_annot.txt')
#don't match b/c some rxns have all 0, such as RXN-4464 which turns 16-HYDROXYPALMITATE into 16-HYDROXYPALMITATE
rxns <- rxns[colnames(s.sp),]
#gapless
s.gapless <- s.sp
s.gapless[s.gapless!=0] <- s.gapless[s.gapless!=0]-10^(-3)
#rename, for below
s.al <- s.gapless
rxns.al <- rxns

##ub
fva.ub <- rep(1000, nrxns); names(fva.ub) <- colnames(s.sp)
#fva.ub0 <- read.csv('s_comb_ub.csv')
#fva.ub[fva.ub0[,1]] <- fva.ub0[,2]

##min meta rxns
#rxns[rxns$nc==1, 'prob'] <- 1.1
rxns.al[unique(c(bm$rxn, vogel$rxn)), 'nc'] <- 1
meta.add.rxns <- read.csv('meta_add_rxns.csv')
jd.rxns <- c(meta.add.rxns$rxn, rxns.al$rxn[rxns.al$nc==1])
rxns.al[jd.rxns, 'prob'] <- 1.1
beta0 <- as.numeric(rxns.al$rxn %in% jd.rxns)

##set list of death conds
#can't yet add sucrose, i think
nut.rxns.lst <- NULL #list('OXYGEN-MOLECULE-TRANS-RXN-L2R')
rad.rxns.ko <- rad.rxns[sapply(rad.rxns, FUN=function(x){ all(x %in% jd.rxns) })]
#set list of rxns
ko.rxns.lst <- c(nut.rxns.lst, rad.rxns.ko)
(n.death.conds <- length(ko.rxns.lst))

##get cutting planes
#cut.met.mat <- cut.mets(s.sp)
cut.met.mat <- Matrix(0, nrow=1, ncol=nrxns)
#cut.rev.mat <- cut.revs(s.sp, rxns)
cut.rev.mat <- Matrix(0, nrow=1, ncol=nrxns)

##instantiate
#conds: vogel & gapless growth, and mult death
#growth cond has nmets+nrxns rows, for sv=0 & v<=M*beta
#decision vars are [v1 mu1 lambda1 mu2 lambda2 ... v2 beta]
#constraints are [vogel.growth death.conds gapless cut.planes obj.constraint gdls]
al.a <- Matrix(0, nrow=nmets+nrxns+2*nrxns*n.death.conds+nmets+nrxns+nrow(cut.met.mat)+nrow(cut.rev.mat)+2, 
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
    rownames(al.a) <- c(paste('fba', c(rownames(s.al), colnames(s.al)), sep='_'), dual.rownames, gapless.rownames, cut.rownames, 'obj', 'gdls')
} else { 
    rownames(al.a) <- c(paste('fba', c(rownames(s.al), colnames(s.al)), sep='_'), gapless.rownames, cut.rownames, 'obj', 'gdls')
}
#colnames
if (n.death.conds>0){ 
    dual.colnames <- paste(rep(c('mu','l'), times=c(nmets, nrxns)), rep(1:n.death.conds, each=nmets+nrxns), '_', c(rownames(s.al), colnames(s.al)), sep='')
    colnames(al.a) <- c(colnames(s.al), dual.colnames, paste('gapless', colnames(s.al), sep='_'), paste('beta', colnames(s.al), sep='_'))
} else { 
    colnames(al.a) <- c(colnames(s.al), paste('gapless', colnames(s.al), sep='_'), paste('beta', colnames(s.al), sep='_')) 
}
#vector names
names(al.lb) <- names(al.ub) <- names(al.obj) <- colnames(al.a)
names(al.rhs) <- names(al.sense) <- rownames(al.a)

##beta bounds & coefficients
al.obj[beta.cols] <- 1.1-rxns.al$prob  
#by default beta lb already 0
al.ub[beta.cols] <- 1
#set betas of desired fluxes to 1
#keep.rxns <- c('biomass', vogel$rxn, unlist(ko.rxns.lst))
#al.lb[paste('beta', unique(keep.rxns), sep='_')] <- 1

##growth on vogel
#sv>=0
#by default, sv rhs already 0
al.a[1:nmets, 1:nrxns] <- s.al
#make 'G' for sv>=0
al.sense[1:nmets] <- 'G'
#Iv-fva.ub*beta<=0
al.a[cbind(nmets + 1:nrxns, 1:nrxns)] <- rep(1, ncol(s.al))
al.a[cbind(nmets + 1:nrxns, beta.cols)] <- -fva.ub
#accidentally tried Iv=1000*beta, w/ good results
al.sense[nmets + 1:nrxns] <- 'L'
#bounds on biomass, fva.ub['biomass']=0.2
al.lb['biomass'] <- 100

##run fba, for use in dual
fba.lp <- FBA(a=s.al, fba.ub=fva.ub)
(fba.obj.val <- fba.lp$obj)
#dual ub: check w/ paschalidis
dual.ub <- rep(1000, nrxns) #fba.obj.val/fva.ub
names(dual.ub)=colnames(s.sp)

##duals
fba.obj <- numeric(nrxns); names(fba.obj) <- colnames(s.sp)
fba.obj['biomass'] <- 1
#for faster assignment
ts.al.nz.ind <- which(t(s.al)!=0)
ts.al.nz <- t(s.al)[ts.al.nz.ind]
#much faster to define in loop instead of defining in fba.dual & passing
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
    al.a[cbind(dci.eq.rows, dci.l.cols)] <- rep(-1, ncol(s.al))
    al.a[cbind(dci.ko.rows, dci.l.cols)] <- rep(1, nrxns)
    
    ##assign betas
    dci.beta.coeff.v <- dual.ub
    #don't want beta=1 -> lambda<0 for ko.rxns or nutrient rxns
    rxns.no.beta <- c(vogel$rxn, ko.rxns)
    dci.beta.coeff.v[which(colnames(s.al) %in% rxns.no.beta)] <- 0
    al.a[cbind(dci.ko.rows, beta.cols)] <- dci.beta.coeff.v
    
    ##bounds on duals
    #ub=Inf by default, above
    al.lb[dci.l.cols] <- -Inf
    #mu>=0 for sv>=0
    al.lb[dci.mu.cols] <- 0
            
    ##bounds & penalties on nut rxns
    al.lb[paste('l', death.cond.i, '_', vogel$rxn, sep='')] <- -Inf
    al.obj[paste('l', death.cond.i, '_', vogel$rxn, sep='')] <- 1000
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
al.a[cbind(gl.ub.rows, gl.v.cols)] <- rep(1, ncol(s.al))
al.a[cbind(gl.ub.rows, beta.cols)] <- rep(-10^4, nrxns)
#rhs
#need to make all nuts highly available to capture feeder pathways for nuts not in vogel
#don't currently have rxns bringing these nuts in extracellular compartment, so allow these extracellular mets to deplete
al.rhs[gl.produce.rows][grep('[e]', names(al.rhs[gl.produce.rows]), fixed=TRUE)] <- -1000
#to get rid of gaplessness constraints for blp, make this -Inf
al.rhs[gl.produce.rows] <- -Inf
#sense
al.sense[gl.produce.rows] <- 'G'
al.sense[gl.ub.rows] <- 'L'

##cutting plane matrices
al.a[nmets+nrxns+2*nrxns*n.death.conds+nmets+nrxns + 1:nrow(cut.met.mat), beta.cols] <- cut.met.mat
al.a[nmets+nrxns+2*nrxns*n.death.conds+nmets+nrxns+nrow(cut.met.mat) + 1:nrow(cut.rev.mat), beta.cols] <- cut.rev.mat
al.sense[nmets+nrxns+2*nrxns*n.death.conds+nmets+nrxns + 1:(nrow(cut.met.mat)+nrow(cut.rev.mat))] <- 'L'

##obj constraint
#this is a minimization
al.a['obj',] <- al.obj
al.sense['obj'] <- 'L'
#start w/ Inf, until have obj value
al.rhs['obj'] <- Inf
al.rhs['obj'] <- al.obj %*% farm.blp$xopt #15245.53

##gdls
beta0 <- numeric(beta.cols)
beta0 <- al.beta>10^-6
al.a['gdls', beta.cols][beta0==1] <- -1
al.a['gdls', beta.cols][beta0==0] <- 1
al.sense['gdls'] <- 'L'
#k-sum(beta0), or Inf to relax
al.rhs['gdls'] <- Inf
al.rhs['gdls'] <- 1-sum(beta0)

##alarm lp
#min weighted penalty of beta's
alarm.lp <- Rcplex(cvec=al.obj, Amat=al.a, bvec=al.rhs, ub=al.ub, lb=al.lb, sense=al.sense)
alarm.lp$stat
x <- alarm.lp$xopt; names(x) <- names(al.obj)

##farm blp
vtype <- rep('C', ncol(al.a))
vtype[beta.cols] <- 'B'
farm.blp <- Rcplex(cvec=al.obj, Amat=al.a, bvec=al.rhs, ub=al.ub, lb=al.lb, sense=al.sense, vtype=vtype, control=list(mipemphasis=0, tilim=11000))
farm.blp$stat
x <- farm.blp$xopt; names(x) <- names(al.obj)
x[beta.cols][!(x[beta.cols] %in% c(0,1))]

##get beta
eps <- 10^-6; thresh <- 10^-6
al.beta <- x[beta.cols]
al.v1 <- x[1:nrxns]
mean(al.beta>1-eps|al.beta<eps)
sum(al.beta>eps)
rxns1 <- sub('beta_', '', names(al.beta)[al.beta>thresh])
beta0 <- as.numeric(al.beta>=thresh)

##vogel growth
fba.r1.ub <- fva.ub[rxns1]
vv <- FBA(s.al[,rxns1], sense='G', fba.ub=fva.ub[rxns1])
#vv <- FBA(s.al[,rxns1], sense='G', goal.rxns=unique(bm$rxn))
#no vogel
#fba.r1.ub["OXYGEN-MOLECULE-TRANS-RXN-L2R"] <- 0
fba.r1.ub[vogel$rxn] <- 0
vn <- FBA(a=s.al[,rxns1], fba.ub=fba.r1.ub, sense='G')
v <- vn$xopt; names(v) <- rxns1

##check dual for a cond
dcond <- 1
for (dcond in 1:length(ko.rxns.lst)){ 
    lam <- x[paste('l', dcond, '_', rxns$rxn, sep='')]
    print(all(lam[paste('l', dcond, '_', vogel$rxn, sep='')]<=10^-6))
}
sum(lam>0 & al.beta==1); lam[lam>0 & al.beta==1]
fd <- fba.dual(a=s.al[,rxns1], ko.rxns=ko.rxns.lst[[dcond]])
fp <- FBA(a=s.al[,rxns1], ko.rxns=ko.rxns.lst[[dcond]])

##compare dual to fba.dual
dcond <- 1
fd <- fba.dual(a=s.al, ko.rxns=union(setdiff(colnames(s.sp), rxns1), ko.rxns.lst[[dcond]]))
#dual
dual.cols <- c(dci.cols, beta.cols)
d.beta.cols <- nmets + nrxns + 1:nrxns
d.obj <- al.obj[dual.cols]; d.obj[d.beta.cols] <- 0
d.lb <- al.lb[dual.cols]; d.ub <- al.ub[dual.cols]
d.lb[d.beta.cols] <- d.ub[d.beta.cols] <- 1
d.a <- al.a[dci.rows,dual.cols]; d.a[1:nrxns, d.beta.cols][cbind(which(beta0==0),which(beta0==0))] <- 0
names(d.obj) <- names(d.ub) <- names(d.lb) <- colnames(al.a)[dual.cols]
fd2 = Rcplex(cvec=d.obj, Amat=d.a, bvec=al.rhs[dci.rows], lb=d.lb, ub=d.ub, sense=al.sense[dci.rows])
cat('Status:', fd2$status, 'w/ dual.obj:', fd2$obj, '\n')
y <- fd2$xopt; names(y) <- colnames(al.a)[dual.cols]; y[paste('l', dcond, '_', rxns.no.beta, sep='')]

##check ko
#rad
ck <- check.ko(s.test=s.al[,rxns1], rxns.test=rxns[rxns1,], ko.lst=rad.rxns, ub=fva.ub[rxns1])
table(ck$obs, ck$pred>eps)
#broad
ck <- check.ko(s.test=s.al[,rxns1], rxns.test=rxns[rxns1,], ko.lst=broad.rxns, obs=1, ub=fva.ub[rxns1])
table(ck$obs, ck$pred>eps)
##write
rr.out <- data.frame(frame=rxns.al[rxns1, 3], flux=x[rxns1], rxns.al[rxns1, -3])
write.table(rr.out, 'rxns_june?.txt', quote=FALSE, row.names=FALSE, sep='\t')

###validate##############################################################################################################
##check feasibility, if stat=5
all(x<=al.ub, x>=al.lb)
sgn <- c('>=', '==', '<=')[match(al.sense, c('G','E','L'))]
check.feas <- apply(cbind(as.numeric(al.a %*% x), sgn, al.rhs), 1, FUN=function(y){ eval(parse(text=paste(y[1], y[2], y[3]))) })

##dual in loop
col.inds <- c(dci.cols, beta.cols)
dual.lp <- Rcplex(cvec=c(al.obj[dci.cols], rep(0, nrxns)), Amat=al.a[dci.rows, col.inds], 
bvec=dual.rhs, sense=dual.sense, lb=c(al.lb[dci.cols], rep(1, length(beta.cols))), ub=al.ub[col.inds])
dual.lp$stat
dual.lp$obj
