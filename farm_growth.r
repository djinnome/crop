##jmd
##8.1.11
##farm_growth.r
##apply farm to incorrectly predicted growth conds, uses sv>=0
#use goal programming, v_bm <= 1

#chmod +x farm_script.r
#call for hour w/ email: bsub -q hour -W 239 -P nc10 -R 'rusage[mem=3]' -N Rscript /msc/neurospora/FBA/farm/farm_growth.r
#priority: bsub -q priority -P nc10 -R 'rusage[mem=3]' -N Rscript /msc/neurospora/FBA/farm/farm_growth.r

##TO DO:
#1. include exports
#2. take only growth conds that we get wrong & grow on universal matrix w/ sv>=0

source('/msc/neurospora/FBA/farm/globals.r')
source('/msc/neurospora/FBA/farm/phenos2rxns.r')
setwd('/msc/neurospora/FBA/farm_data')
options(stringsAsFactors=FALSE)

##read
filt.rxns <- read.delim('Neurospora/nc10.filtered')
wrong.dir.rxns <- read.delim('nc_wrong_dir_rxns.txt')
s.sp <- read.S('s_bal.txt')
nmets <- nrow(s.sp); nrxns <- ncol(s.sp)
#add rxns
meta.add.rxns <- read.delim('meta_add_rxns.txt')
filt.add.rxns <- read.delim('filt_add_rxns.txt')
wrong.dir.add.rxns <- read.delim('wrong_dir_add_rxns.txt')
#rxns annot
rxns <- read.delim('rxn_bal_annot.txt')
#in case don't match
rxns <- rxns[colnames(s.sp),]
print('Check that annotations match S matrix')
all(rownames(rxns)==colnames(s.sp))
#gapless
s.gapless <- s.sp
s.gapless[s.gapless!=0] <- s.gapless[s.gapless!=0]-10^(-3)
#rename, for below
s.al <- s.gapless
rxns.al <- rxns
#ub
fva.ub <- rep(1000, nrxns); names(fva.ub) <- colnames(s.al)

##set rxns
rid.rxns <- c('SUCROSE-SYNTHASE-RXN-CPD-12575/BETA-D-FRUCTOSE//SUCROSE/UDP/PROTON.46.-R2L', 'XANPRIBOSYLTRAN-RXN-R2L',
'FRUCTOSE-6-PHOSPHATE-PHOSPHOKETOLASE-RXN-L2R', 'GLUCONOKIN-RXN-L2R', "ARYLFORMAMIDASE-RXN-L2R")
jd.rxns <- setdiff(union(c(rownames(meta.add.rxns$rxn), rownames(filt.add.rxns)), setdiff(rownames(rxns.al)[rxns.al$nc==1], filt.rxns$rxn)), rid.rxns)
jd.rxns <- jd.rxns[!is.na(jd.rxns)]

##check ko
ck0 <- check.ko(s.test=s.al[,jd.rxns], rxns.test=rxns.al[jd.rxns,], ko.lst=broad.rxns, ub=fva.ub[jd.rxns], print.v=FALSE)
table(ck0$obs, ck0$pred>10^-6)
#write.table(ck0, 'check_ko_broad_aug1.txt', sep='\t', quote=FALSE, row.names=FALSE)

##set list of conds
ko.rxns.lst <- broad.rxns[!is.na(ck0$pred) & ck0$pred<10^-6]
(n.conds <- length(ko.rxns.lst))

##instantiate
#decision vars are [v1 v2 v3 ... beta]
#add row for gdls
al.a <- Matrix(0, nrow=(nmets+nrxns)*n.conds+1, ncol=nrxns*(n.conds+1))
al.lb <- al.obj <- rep(0, ncol(al.a))
al.ub <- c(rep(1000, n.conds*nrxns), rep(1, nrxns))
al.rhs <- rep(0, nrow(al.a))
al.sense <- rep('G', nrow(al.a))
#index beta
beta.cols <- nrxns*n.conds + 1:nrxns

##names
#rownames
rownames(al.a) <- c(paste(rownames(s.sp), '_cond', rep(1:n.conds, each=nrow(s.sp)), sep=''), 
paste(colnames(s.sp), '_cond', rep(1:n.conds, each=ncol(s.sp)), sep=''), 'gdls')
#colnames
colnames(al.a) <- c(paste(colnames(s.sp), '_cond', rep(1:n.conds, each=ncol(s.sp)), sep=''), 
paste('beta_', colnames(s.sp), sep=''))
#vector names
names(al.lb) <- names(al.ub) <- names(al.obj) <- colnames(al.a)
names(al.rhs) <- names(al.sense) <- rownames(al.a)

##beta bounds & coefficients
rxns.al$prob[rxns.al$prob<0] <- -0.05
rxns.al[jd.rxns, 'prob'] <- 1.1
#set obj coeffs: set to not get rid of any, since want min change and don't have gaplessness
al.obj[beta.cols] <- 1.1-rxns.al$prob
#by default beta lb already 0
al.ub[beta.cols] <- 1
#goal prog of bm
al.ub[paste('biomass_cond', 1:n.conds, sep='')] <- 1
al.obj[paste('biomass_cond', 1:n.conds, sep='')] <- -10

##fba
#format of al.a is [s 0; 0 s 0;...; -I 0 1000*I; 0 -I 0 1000*i]
for (cond.ind in 1:n.conds){
    ci.cols <- nrxns*(cond.ind-1) + 1:nrxns
    #sv>=0
    al.a[nmets*(cond.ind-1) + 1:nmets, ci.cols] <- s.al
    #-Iv+1000*beta>=0
    ci.vlim.rows <- nmets*n.conds+nrxns*(cond.ind-1) + 1:nrxns
    al.a[ci.vlim.rows, ci.cols] <- -1*Diagonal(ncol(s.sp))
    al.a[ci.vlim.rows, beta.cols] <- Diagonal(x=fva.ub)
    #ko
    al.ub[paste(ko.rxns.lst[[cond.ind]], '_cond', cond.ind, sep='')] <- 0    
    #print
    print(paste('added condition', cond.ind))
}

###lp####################################################################################################################
##alarm
alarm.lp <- Rcplex(cvec=al.obj, Amat=al.a, bvec=al.rhs, ub=al.ub, lb=al.lb, sense=al.sense, objsense='min')
alarm.lp$stat
x <- alarm.lp$x; names(x) <- colnames(al.a)
x[paste('biomass_cond', 1:n.conds, sep='')]

##get beta
eps <- 10^-6; thresh <- 10^-9
al.beta <- x[beta.cols]; names(al.beta) <- colnames(s.al)
(new.rxns <- setdiff(names(al.beta)[al.beta>thresh], jd.rxns))
rxns1 <- union(names(al.beta)[al.beta>thresh], jd.rxns)
beta0 <- as.numeric(names(al.beta) %in% rxns1)

##check ko
#broad
ckb <- check.ko(s.test=s.al[,rxns1], rxns.test=rxns.al[rxns1,], ko.lst=ko.rxns.lst, obs=1, ub=fva.ub[rxns1])
table(ckb$obs, ckb$pred>eps)
#rad
ckr <- check.ko(s.test=s.al[,rxns1], rxns.test=rxns.al[rxns1,], ko.lst=rad.rxns, obs=1, ub=fva.ub[rxns1])
table(ckr$obs, ckr$pred>eps)

###blp###################################################################################################################
##gdls
al.a['gdls', beta.cols][beta0==1] <- -1
al.a['gdls', beta.cols][beta0==0] <- 1
al.sense['gdls'] <- 'L'
#k-sum(beta0), or Inf to relax
n.change <- 10
al.rhs['gdls'] <- n.change-sum(beta0)
print(paste('N change =', n.change))

##farm
vtype <- rep('C', ncol(al.a)); vtype[beta.cols] <- 'B'
farm.blp <- Rcplex(cvec=al.obj, Amat=al.a, bvec=al.rhs, ub=al.ub, lb=al.lb, sense=al.sense, vtype=vtype, control=list(mipemphasis=0, tilim=1*3600))
farm.blp$stat
x <- farm.blp$xopt; names(x) <- names(al.obj)
x[paste('biomass_cond', 1:n.conds, sep='')]
print("Non-binary CPLEX 'binary' variables")
#y <- x[beta.cols][!(x[beta.cols] %in% c(0,1))]; y[order(-y)]

##get beta
eps <- 10^-6; thresh <- 10^-9
al.beta <- x[beta.cols]; names(al.beta) <- colnames(s.al)
(new.rxns <- setdiff(names(al.beta)[al.beta>thresh], jd.rxns))
#new0 <- c("H2NTPEPIM-RXN-L2R", "RXN-10857-L2R", "RXN-1725-L2R", "RXN0-6368-L2R")
rxns1 <- union(names(al.beta)[al.beta>thresh], jd.rxns)
beta0 <- as.numeric(names(al.beta) %in% rxns1)

##check ko
ckb <- check.ko(s.test=s.al[,rxns1], rxns.test=rxns.al[rxns1,], ko.lst=broad.rxns, obs=1, ub=fva.ub[rxns1])
table(ckb$obs, ckb$pred>eps)
