##jmd
##6.23.11
##min_meta_rxns.r

source('/msc/neurospora/FBA/farm/farm_config.r')

##make rxns
rxns <- data.frame(matrix(0, nrow=ncol(s.sp), ncol=2, dimnames=list(colnames(s.sp), c('nc', 'prob'))))
rxns[rownames(rxns) %in% c(nc$rxn, bm$rxn, bm.goals$rxn, vogel$rxn), 'nc'] <- 1
#check
print('Check that annotations match matrix')
all(colnames(s.sp)==rownames(rxns))

##rename, for below
s.al <- sg
rxns.al <- rxns

##ub
fva.ub <- rep(1000, nrxns); names(fva.ub) <- colnames(s.al)

##check nc growth
#goal.lp <- FBA(a=s.al[,rxns.al$nc==1], fba.obj.rxns=NULL, goal.rxns=unique(bm.goals$rxn))
#v <- goal.lp$xopt; names(v) <- colnames(s.al)[rxns.al$nc==1]
#vg <- v[unique(bm.goals$rxn)]; vg[vg==0]

##min bad rxns
rxns.al[rxns.al$nc==1, 'prob'] <- 1.1
#filtered rxns
#sum(filt.rxns$rxn %in% rownames(rxns.al))
rxns.al[rownames(rxns.al) %in% filt.rxns$rxn, 'prob'] <- 0.9
#wrong dir rxns
rxns.al[rownames(rxns.al) %in% wrong.dir.rxns$rxn, 'prob'] <- 1

##imbalanced rxns
bal.mat <- read.csv('rxns_imbalance.csv', row.names=1)
bal.mat <- bal.mat[,c('C', 'O', 'N', 'P', 'S')]
unb.rxns <- rownames(bal.mat)[rowSums(bal.mat>0)>0]
rxns.al[rownames(rxns.al) %in% unb.rxns & rxns.al$nc==1, 'prob'] <- -50
rxns.al[rownames(rxns.al) %in% unb.rxns & rxns.al$nc==0, 'prob'] <- -50
#very naughty rxn
rxns.al['NAD+-ADP-RIBOSYLTRANSFERASE-RXN-R2L', 'prob'] <- -10**3

##decision vars are [v beta]
#a = [s 0; I -1000]
obj <- c(numeric(nrxns), 1.1-rxns.al$prob)
amat <- rBind(cBind(s.al, Matrix(0, nrow=nrow(s.al), ncol=nrxns)), cBind(Diagonal(n=nrxns), -1000*Diagonal(n=nrxns)))
rhs <- numeric(nmets+nrxns)
ub <- rep(c(10**3, 1), times=c(nrxns, nrxns))
lb <- numeric(2*nrxns)
sense <- rep(c('G', 'L'), times=c(nmets, nrxns))
#names
names(lb) <- names(ub) <- names(obj) <- c(colnames(s.al), paste('beta_', colnames(s.al), sep=''))
lb['biomass'] <- 1
#lp
min.meta.lp <- Rcplex(cvec=obj, Amat=amat, bvec=rhs, ub=ub, lb=lb, sense=sense, control=list(trace=0))
print('CPLEX status of LP')
min.meta.lp$stat
x <- min.meta.lp$xopt; names(x) <- names(obj)
v <- x[1:nrxns]
beta <- x[nrxns + 1:nrxns]

##get min adjust
#threshold for beta of meta rxns
for (thresh in c(10**-6, 10**-9, 10**-15, 0)){
    rr <- rxns.al[beta>thresh,]
    #check goals of beta>thresh
    #fva.ub[beta<thresh] <- 0
    #goals <- FBA(a=s.al, fba.obj.rxns=NULL, goal.rxns=unique(bm.goals$rxn), fba.ub=fva.ub)
    #get needed meta rxns, nrows are number of rxns
    cat('Number of needed meta rxns:', nrow(meta.add.rxns <- rr[rr$nc==0,]), '\n')
    cat('Number of needed filtered rxns:', nrow(filt.add.rxns <- rr[rownames(rr) %in% filt.rxns$rxn,]), '\n')
    cat('Number of needed wrong direction rxns:', nrow(wrong.dir.add.rxns <- rr[rownames(rr) %in% wrong.dir.rxns$rxn,]), '\n')
    #check growth of nc & meta.rxns combined
    #rxns1 <- c(rownames(meta.rxns), rownames(rxns.al)[rxns.al$nc==1])
    rxns1 <- union(rownames(rxns.al)[beta>thresh], setdiff(rownames(rxns.al)[rxns.al$nc==1], union(wrong.dir.rxns$rxn, filt.rxns$rxn)))
    print(paste('Testing if new rxn set gives FBA growth w/ threshold', thresh))
    fba <- FBA(s.al[,rxns1], control=list(trace=0))
    cat('Unbalanced rxns to be added:', unb.add <- intersect(unb.rxns, c(rownames(meta.add.rxns), rownames(filt.add.rxns), rownames(wrong.dir.add.rxns))), '\n')
    if (fba$obj>10**-6){
        fba.v <- fba$xopt; names(fba.v) <- rxns1
        flux.mat <- print.flux(v=fba.v, out.file='fluxes.tsv')
        #write these reaction
        cat('writing files \n')
        write.table(meta.add.rxns, 'meta_add_rxns.txt', quote=FALSE, sep='\t')
        write.table(filt.add.rxns, 'filt_add_rxns.txt', quote=FALSE, sep='\t')
        write.table(wrong.dir.add.rxns, 'wrong_dir_add_rxns.txt', quote=FALSE, sep='\t')
        break #out of the for loop
    }#end if
}#end for

#(fpr <- rr[rownames(rr) %in% free.prod.rxns,])

##check ko
#source('/msc/neurospora/FBA/farm/phenos2rxns.r')
#rad
#ck <- check.ko(s.test=s.al[,rxns1], rxns.test=rxns.al[rxns1,], ko.lst=rad.rxns, ub=fva.ub[rxns1]); table(ck$obs, ck$pred>10^-6)
#broad
#ckb <- check.ko(s.test=s.al[,rxns1], rxns.test=rxns.al[rxns1,], ko.lst=broad.rxns, obs=1, ub=fva.ub[rxns1]); table(ckb$obs, ckb$pred>10^-6)
