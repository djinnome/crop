##jmd
##7.25.11
##farm_script.r
##apply farm to meta db

#chmod +x farm_script.r
#call for hour w/ email: bsub -q hour -W 239 -P nc10 -R 'rusage[mem=3]' -N Rscript /msc/neurospora/FBA/farm/farm_script.r
#priority: bsub -q priority -P nc10 -R 'rusage[mem=3]' -N Rscript /msc/neurospora/FBA/farm/farm_script.r

source('/msc/neurospora/FBA/farm/farm_config.r')
source('/msc/neurospora/FBA/farm/phenos2rxns.r')

##read
#rxns annot
rxns <- read.delim('rxn_annot.txt')
rxns <- rxns[colnames(s.sp),]
all(rownames(rxns)==colnames(s.sp))

##subset for radford ko's
model.rxns <- read.csv('model_rxns.csv')[,1]
ss <- intersect(model.rxns, colnames(sg))
s.al <- sg[,ss]
s.al <- s.al[rowSums(abs(s.al))>0,]
rxns.al <- rxns[ss,]
nmets <- nrow(s.al); nrxns <- ncol(s.al)

##get rid of no flux rxns
#prune dead-end mets in reversible model which are connected to rest of network by a single reversible path, so can carry flux in irrev model
s.rev <- t(apply(s.al, 1, FUN=function(x){
    tapply(x, INDEX=gsub('-L2R|-R2L', '', colnames(s.al)), FUN=function(y){ sum(abs(y)) })
}))
(deadend.mets <- setdiff(rownames(s.rev)[rowSums(s.rev)<=1], grep("CCO-EXTRACELLULAR", rownames(s.al), value=TRUE)) )

#max.v <- max.flux(a=s.al, v.ub=1)
max.nv <- max.n.flux(s.al, se='E', allow.all.trans=TRUE)
#find no flux rxns under sv=0, st can't cheat thru "export dilution"
max.nv2 <- max.n.flux(s.sp[rownames(s.al),colnames(s.al)], se='E', allow.all.trans=TRUE)
no.v <- names(max.nv)[max.nv<=0]
no.v2 <- names(max.nv2)[max.nv2<=0]
#write.table(rxns.al[no.v, c('FRAME.ID', 'rxn')], 'no_flux_rxns.txt', quote=FALSE, row.names=FALSE)
#no.v <- read.table('no_flux_rxns.txt', header=TRUE)
s.al <- s.al[,!(colnames(s.al) %in% no.v)]
s.al <- s.al[rowSums(abs(s.al))>0,]
rxns.al <- rxns.al[colnames(s.al),]
#summarize
(nmets <- nrow(s.al)); (nrxns <- ncol(s.al))
beta0 <- rep(1, nrxns)

##ub
#DO *NOT* change this ub to Inf, it is used in al.a
fva.ub <- rep(10**3, nrxns); names(fva.ub) <- colnames(s.al)
fva.ub[intersect(vogel$rxn, colnames(s.al))] <- 1000; fva.ub['SUCROSE-TRANS-RXN-L2R'] <- 5
#fba
cat('Checking growth:', '\n')
f <- FBA(s.al, sense='E', fba.ub=fva.ub)
f <- FBA(s.al, sense='E', fba.ub=fva.ub, fba.obj='FullBiomassComposition')
#trans.rxns doesn't include vogels, since their rxn names start w/ cpd name
#trans.rxns <- grep('^TRANS-RXN', colnames(s.al), value=TRUE)
#fva.ub[trans.rxns] <- 10^6
#fva.ub0 <- read.csv('s_comb_ub.csv')
#fva.ub[fva.ub0[,1]] <- fva.ub0[,2]

##set list of death conds
#only use death rxns that are not correctly predicted
rad.rxns.ko <- rad.rxns[!is.na(ck0$pred) & ck0$pred>10^-6]
#set list of rxns
#ko.rxns.lst <- rad.rxns.ko
(ko.rxns.lst <- rad.rxns['NCU09111.5'])
(n.death.conds <- length(ko.rxns.lst))

##get cutting planes
#these only allow rxns which can carry flux on vogels
cut.met.mat <- cut.mets(s.al, sense='E')
cut.rev.mat <- cut.revs(s.al, rxns.al)
#cut.met.mat <- Matrix(0, nrow=1, ncol=nrxns)
#cut.rev.mat <- Matrix(0, nrow=1, ncol=nrxns)

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
#set rnxs of value for dual
#c.rxns <- c('SUCROSE-TRANS-RXN-L2R','CIT-TRANS-RXN-L2R')
c.rxns <- 'biomass'

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
#set coeff for unbalanced rxns
rxns.al$prob[rxns.al$prob<0] <- -0.05
#set obj coeffs: set to not get rid of any, since want min change and don't have gaplessness
al.obj[beta.cols] <- -0.1-rxns.al$prob
#by default beta lb already 0
al.ub[beta.cols] <- 1

##set betas of desired fluxes to 1
#already got rid of no-flux rxns, making this faster
#perform single min |v| to further reduce num of rxns
mv.lb <- numeric(nrxns); names(mv.lb) <- colnames(s.al); mv.lb['biomass'] <- 1
min.abs.v <- Rcplex(cvec=rep(1, nrxns), Amat=s.al, bvec=numeric(nmets), sense='E', lb=mv.lb)
nec.rxns0 <- colnames(s.al)[min.abs.v$xopt>0]
#min.v <- min.flux(a=s.al, rxns.min=nec.rxns0, min.bm=1)
#need.rxns <- names(min.v[min.v>0])
#write.table(need.rxns, 'necessary_rxns.txt', quote=FALSE, row.names=FALSE)
need.rxns <- read.table('necessary_rxns.txt', header=TRUE)[,1]
keep.rxns <- c('biomass', vogel$rxn, grep('^TRANS-RXN2T-232-', colnames(s.al), value=TRUE), need.rxns, "RXN0-5114-L2R", "TRANS-RXN2T-81-L2R[CCO-PM-FUNGI]")
al.lb[intersect(paste('beta', unique(keep.rxns), sep='_'), names(al.lb))] <- 1

##growth on vogel
#by default, sv rhs already 0
al.a[1:nmets, 1:nrxns] <- s.al
#make 'E' for sv=0
al.sense[1:nmets] <- 'E'
#Iv-fva.ub*beta<=0
al.a[cbind(nmets + 1:nrxns, 1:nrxns)] <- rep(1, ncol(s.al))
al.a[cbind(nmets + 1:nrxns, beta.cols)] <- -fva.ub
#accidentally tried Iv=1000*beta, w/ good results
al.sense[nmets + 1:nrxns] <- 'L'
#bounds on biomass, fva.ub['biomass']=0.2
al.lb['biomass'] <- 0.1

##get dual ub's
#use dual.ub=1 when only limited flux is biomass
dual.ub <- rep(1, nrxns)
names(dual.ub) <- colnames(s.al)
#need.rxns won't be KO'd, so their lambda's<=0
dual.ub[intersect(keep.rxns, names(dual.ub))] <- 0

##duals
fba.obj <- numeric(nrxns); names(fba.obj) <- colnames(s.al); fba.obj['biomass'] <- 1
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
    rxns.no.beta <- c(c.rxns, ko.rxns)
    dci.beta.coeff.v[which(colnames(s.al) %in% rxns.no.beta)] <- 0
    al.a[cbind(dci.ko.rows, beta.cols)] <- dci.beta.coeff.v
    
    ##bounds on duals
    #ub=Inf by default, above
    al.lb[dci.l.cols] <- -Inf
    #mu free  for sv=0
    al.lb[dci.mu.cols] <- -Inf
            
    ##bounds & penalties on nut rxns
    al.lb[paste('l', death.cond.i, '_', c.rxns, sep='')] <- 0
    #al.obj[paste('l', death.cond.i, '_', vogel$rxn, sep='')] <- 100
    al.obj[paste('l', death.cond.i, '_', c.rxns, sep='')] <- 100
    
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
#al.rhs['obj'] <- 1244 #al.obj %*% farm.blp$xopt

##gdls
al.a['gdls', beta.cols][beta0==1] <- -1
al.a['gdls', beta.cols][beta0==0] <- 1
al.sense['gdls'] <- 'L'
#k-sum(beta0), or Inf to relax
n.ko <- Inf
al.rhs['gdls'] <- n.ko-sum(beta0)
print(paste('N KO =', n.ko))

##farm blp
vtype <- rep('C', ncol(al.a))
vtype[beta.cols] <- 'B'
all(!is.na(al.obj), is.finite(fva.ub))

farm.blp <- Rcplex(cvec=al.obj, Amat=al.a, bvec=al.rhs, ub=al.ub, lb=al.lb, sense=al.sense, vtype=vtype, control=list(mipemphasis=0, tilim=3600/10)) #epagap=100*(n.death.conds-1)+25
farm.blp$stat
x <- farm.blp$xopt; names(x) <- names(al.obj)
print("Non-binary CPLEX 'binary' variables")
y <- x[beta.cols][!is.na(x[beta.cols]) & !(x[beta.cols] %in% c(0,1))]; y[order(-y)]

##putative consistent radford predictions
l.bm <- x[paste('l', 1:n.death.conds, '_biomass', sep='')]
print('Duals predicted to have reduced growth.')
l.bm[l.bm<1-10^-6]
#ko.rxns.lst[l.bm<1-10^-6]
#get >500 lambda's>0 when use dual.ub=1000
#ll <- x[paste('l1_', colnames(s.al), sep='')]; ll <- ll[order(-ll)]; ll[1:9]

##get beta
eps <- 10^-6; thresh <- 10^-6
al.beta <- x[beta.cols]
al.v1 <- x[1:nrxns]
sum(al.beta>eps)
rxns1 <- sub('beta_', '', names(al.beta)[al.beta>thresh])
beta0 <- as.numeric(al.beta>=thresh)

#print(paste("rid rxns:", setdiff(jd.rxns, rxns1)))
(rid2 <- setdiff(colnames(s.al), rxns1))
#print(paste("rid rxns = c(", paste(rid2, collapse="\', \'"), ")", sep="\'"))

##check ko
f <- FBA(s.al[,rxns1], sense='E')
#rad
ck <- check.ko(s.test=s.al[,rxns1], ko.lst=ko.rxns.lst, sense='E', annot.df=rad0[names(rad.rxns),])
summary(as.factor(ck$pred>10^-6))
ck[ck$obs==0 & ck0$pred>0,]
#broad
#ckb <- check.ko(s.test=s.al[,rxns1], rxns.test=rxns.al[rxns1,], ko.lst=broad.rxns, obs=1, ub=fva.ub[rxns1])
#table(ckb$obs, ckb$pred>eps)

##write
#write.csv(model.rxns, 'model_rxns.csv', row.names=FALSE)
#rr.out <- data.frame(frame=rxns.al[,3], flux=al.v1, beta=al.beta, name=rownames(rxns.al))
#write.table(rr.out, paste('rxns_script_aug1_nKO', n.ko, '.txt', sep=''), quote=FALSE, row.names=FALSE, sep='\t')
