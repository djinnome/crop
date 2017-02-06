##jmd
##3.23.11
##alarm_growth.r
##universal matrix w/ growth only in 1 media, list of non-lethal KOs, no imports for non-media, & sv>=0

##to do
#need fva
#throttle primal, and constrain opt >= old opt - (0.9-rho), to hasten lp soln

setwd('/msc/neurospora/FBA/farm')

source('check_ko.r')
source('check_fba.r')

##get experimental KO's

##set list of ko rxns as growth conds
#NULL means no ko
ko.rxns.lst <- list(NULL)
(n.conds <- length(ko.rxns.lst))

##read
rxns <- read.delim('rxn_annot.txt')
bm <- read.delim('biomass_s.txt')
vogel <- read.delim('vogel_trans_s.txt')
s <- read.delim('s15.txt', as.is=FALSE, na='NIL')
s.sp <- sparseMatrix(i=as.numeric(s$compound), j=as.numeric(s$rxn), x=as.real(s$coeff))
dimnames(s.sp) <- list(levels(s$compound), levels(s$rxn))
nmets <- nrow(s.sp); nrxns <- ncol(s.sp)
#ub
fba.ub0 <- read.csv('s_comb_ub.csv')
fba.ub <- fba.ub0[,2]; names(fba.ub) <- fba.ub0[,1]

##instantiate
#decision vars are [v1 v2 v3 ... beta]
al.a <- Matrix(0, nrow=(nmets+nrxns)*n.conds, ncol=nrxns*(n.conds+1))
al.lb <- al.obj <- rep(0, ncol(al.a))
#can't have > 10^3 mmol/g
al.ub <- c(fba.ub, rep(1, nrxns))
al.rhs <- rep(0, nrow(al.a))
al.sense <- rep('G', nrow(al.a))
#index beta
beta.cols <- nrxns*n.conds + 1:nrxns
#set carbon source rnxs
c.rxns <- c('SUCROSE-TRANS-RXN-L2R','CIT-TRANS-RXN-L2R')

##names
#rownames
rownames(al.a) <- c(paste(rownames(s.sp), '_cond', rep(1:n.conds, each=nrow(s.sp)), sep=''), 
paste(colnames(s.sp), '_cond', rep(1:n.conds, each=ncol(s.sp)), sep=''))
#colnames
colnames(al.a) <- c(paste(colnames(s.sp), '_cond', rep(1:n.conds, each=ncol(s.sp)), sep=''), 
paste('beta_', colnames(s.sp), sep=''))
#vector names
names(al.lb) <- names(al.ub) <- names(al.obj) <- colnames(al.a)
names(al.rhs) <- names(al.sense) <- rownames(al.a)

##fba
#format of al.a is [s 0; 0 s 0;...; -I 0 1000*I; 0 -I 0 1000*i]
for (cond.ind in 1:n.conds){
    ci.cols <- nrxns*(cond.ind-1) + 1:nrxns
    #sv>=0
    al.a[nmets*(cond.ind-1) + 1:nmets, ci.cols] <- s.sp
    #-Iv+1000*beta>=0
    ci.vlim.rows <- nmets*n.conds+nrxns*(cond.ind-1) + 1:nrxns
    al.a[ci.vlim.rows, ci.cols] <- -1*Diagonal(ncol(s.sp))
    al.a[ci.vlim.rows, beta.cols] <- Diagonal(x=fba.ub)
    
    #bounds on fluxes
    #if this ub=100, growth ~= 3
    #al.ub[paste(c.rxns, '_cond', cond.ind, sep='')] <- 10
    al.lb[paste('biomass_cond', cond.ind, sep='')] <- 0.1
    #ko
    if (!is.null(ko.rxns.lst[[cond.ind]])){
        al.ub[paste(ko.rxns.lst[[cond.ind]], '_cond', cond.ind, sep='')] <- 0
    }
}

##beta's
probs <- rxns$prob
names(probs) <- rxns$rxn
#probs[c('biomass', vogel$rxn)] <- 1
al.obj[beta.cols] <- probs-0.9
#by default beta lb already 0
#al.ub[beta.cols] <- 1

##alarm lp
#max probs-0.9 of beta's
alarm.lp <- Rcplex(cvec=al.obj, Amat=al.a, bvec=al.rhs, ub=al.ub, lb=al.lb, sense=al.sense, objsense='max')
alarm.lp$stat

##extract soln
x <- alarm.lp$xopt; names(x) <- names(al.obj)
al.beta <- x[beta.cols]
mean(al.beta>0.9|al.beta<0.1)
sum(al.beta>10^(-6))
rxns1 <- sub('beta_', '', names(al.beta)[al.beta>0])

##test
#vogel growth
vv <- fba(s.sp[,rxns1], sense='G')
#check ko
ck <- check.ko(s.sp[,rxns1], rxns[rxns1,])
table(ck$obs, ck$pred)

##write
rr.out <- data.frame(frame=rxns[rxns1, 3], flux=x[rxns1], rxns[rxns1, -3])
#write.table(rr.out, 'rxns_beta1_svg0.txt', quote=FALSE, row.names=FALSE, sep='\t')

###validate##############################################################################################################
