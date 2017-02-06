##jmd
##6.4.11
##farm_vogel_seed.r
##universal matrix w/ growth only in vogel, no KO, no imports for non-media, & sv>=0

##to do
#throttle primal, and constrain opt >= old opt - (0.9-rho), to hasten lp soln

library('Matrix')
library('Rcplex')

setwd('/msc/neurospora/FBA/farm')
source('check_ko.r'); source('check_fba.r'); source('fba.r')
source('cut_mets.r'); source('cut_revs.r')
setwd('/msc/neurospora/FBA/seed')

##read
bm <- read.delim('biomass_seed.txt')
rxns <- read.delim('rxns_seed_annot.txt')
vogel <- read.delim('vogel_seed.txt')
s <- read.delim('s_seed.txt', as.is=FALSE, na='NIL')
s.sp <- sparseMatrix(i=as.numeric(s$compound), j=as.numeric(s$rxn), x=as.real(s$coeff))
dimnames(s.sp) <- list(levels(s$compound), levels(s$rxn))
nmets <- nrow(s.sp); nrxns <- ncol(s.sp)
#ub
fba.ub <- rep(1000, nrxns)
names(fba.ub) <- colnames(s.sp)

##get cutting planes
cut.met.mat <- cut.mets(s.sp)
#cut.met.mat <- Matrix(0, nrow=1, ncol=nrxns)
#cut.rev.mat <- cut.revs(s.sp, rxns)
cut.rev.mat <- Matrix(0, nrow=1, ncol=nrxns)

##instantiate
#decision vars are [v beta]
al.a <- Matrix(0, nrow=nmets+nrxns+nrow(cut.met.mat)+nrow(cut.rev.mat), ncol=2*nrxns)
al.lb <- al.obj <- rep(0, ncol(al.a))
#fba.ub is coefficient of beta in matrix
al.ub <- rep(c(1000, 1), each=nrxns)
al.rhs <- rep(0, nrow(al.a))
al.sense <- c(rep('G', nmets), rep('L', nrow(al.a)-nmets))
#index beta
beta.cols <- nrxns + 1:nrxns
#set carbon source rnxs
c.rxns <- c('SUCROSE-TRANS-RXN-L2R','CIT-TRANS-RXN-L2R')

##names
cut.rownames <- c(paste('cut_met', 1:nrow(cut.met.mat), sep=''), paste('cut_rev', 1:nrow(cut.rev.mat), sep=''))
#rownames
rownames(al.a) <- c(rownames(s.sp), colnames(s.sp), cut.rownames)
#colnames
colnames(al.a) <- c(colnames(s.sp), paste('beta_', colnames(s.sp), sep=''))
#vector names
names(al.lb) <- names(al.ub) <- names(al.obj) <- colnames(al.a)
names(al.rhs) <- names(al.sense) <- rownames(al.a)

##growth on vogel
#format of al.a is [s 0; I -M*I; 0 cut.planes]
v.cols <- 1:nrxns
#sv>=0
al.a[1:nmets, v.cols] <- s.sp
#Iv-1000*beta<=0
vlim.rows <- nmets+1:nrxns
al.a[vlim.rows, v.cols] <- Diagonal(nrxns)
al.a[vlim.rows, beta.cols] <- -1*Diagonal(x=fba.ub)

##cutting plane matrices
al.a[nmets+nrxns + 1:nrow(cut.met.mat), beta.cols] <- cut.met.mat
al.a[nmets+nrxns+nrow(cut.met.mat)+ 1:nrow(cut.rev.mat), beta.cols] <- cut.rev.mat

##beta's
probs <- rxns$prob
names(probs) <- rownames(rxns)
al.obj[beta.cols] <- probs-0.99

#bounds on fluxes
#if this ub=100, growth ~= 3
#al.ub[paste(c.rxns, '_cond', cond.ind, sep='')] <- 10
al.lb['biomass'] <- 1
#on some beta's
al.lb[paste('beta_', c('biomass', c.rxns), sep='')] <- 1

##alarm lp
#max probs-0.9 of beta's
alarm.lp <- Rcplex(cvec=al.obj, Amat=al.a, bvec=al.rhs, ub=al.ub, lb=al.lb, sense=al.sense, objsense='max')
alarm.lp$stat

##extract soln
eps <- 10^-6
x <- alarm.lp$xopt; names(x) <- names(al.obj)
v <- x[1:nrxns]
al.beta <- x[beta.cols]; names(al.beta) <- sub('beta_', '', names(al.beta))
mean(al.beta>1-eps|al.beta<eps)
sum(al.beta>eps)
rxns1 <- names(al.beta)[al.beta>0]

##test
#vogel growth
ss <- s.sp[,rxns1]
ss[ss!=0] <- ss[ss!=0]-10^-6
vv <- fba(ss, sense='G', fba.ub=fba.ub[rxns1])

##check ko
ko.rad <- read.delim('/msc/neurospora/FBA/Neurospora/eCompendium/phenotype-experiment.txt', as.is=TRUE)
#get rid of 3-level ec
ko.rad <- ko.rad[ko.rad$Nut!='P-AMINO-BENZOATE',]
ko.rad <- ko.rad[!duplicated(ko.rad$EC),]
#run
ck <- check.ko.ec(s.sp[,rxns1], rxns[rxns1,], ko=ko.rad)
table(ck$obs, ck$pred>eps)

##throttle
length(thr <- which(al.beta<10^-1 & al.beta>eps))
#fba.ub2 <- fba.ub
#fba.ub2[thr] <- fba.ub2[thr]*al.beta[thr]
#al.a[vlim.rows, beta.cols] <- -1*Diagonal(x=fba.ub2)

##write
rr.out <- data.frame(frame=rxns[rxns1, 3], flux=x[rxns1], rxns[rxns1, -3])
write.table(rr.out, 'rxns_may25.txt', quote=FALSE, row.names=FALSE, sep='\t')

###validate##############################################################################################################
probs[c('biomass', vogel$rxn)] <- 1
#by default beta lb already 0
al.ub[beta.cols] <- 1
