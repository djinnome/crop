##jmd
##8.16.11
##run_check_ko.r

source('/msc/neurospora/FBA/farm/farm_config.r')
source('/msc/neurospora/FBA/farm/phenos2rxns.r')

##subset
rck.rxns <- read.csv('model_rxns.csv')[,1]
s.ck <- ss.mat(sg, c=intersect(rck.rxns, colnames(sg)))
sp.ck <- s.sp[rownames(s.ck), colnames(s.ck)]
#ub
fva.ub <- rep(10**3, ncol(s.ck)); names(fva.ub) <- colnames(s.ck); fva.ub['SUCROSE-TRANS-RXN-L2R'] <- 5

##wt growth
cat('\n Wild-type growth on vogels: \n')
fba.lp <- FBA(s.ck, sense='E', fba.ub=fva.ub, control=list(trace=0)); v <- fba.lp$xopt
cat('\n Wild-type FullBiomassComposition growth on vogels: \n')
fba.full.lp <- FBA(s.ck, sense='E', fba.ub=fva.ub, fba.obj.rxn='FullBiomassComposition', control=list(trace=0))

##KO msg
cat('\n Checking KOs: TRUE=growth; FALSE=no growth\n')

##broad
ckb <- call.check.ko(s.ck, broad.rxns, set.name='Broad', obs=1, cprim=numeric(ncol(s.ck)))

##check radford ko
#nut.rxns <- list(o2='OXYGEN-MOLECULE-TRANS-RXN-L2R', c.rxns=c('CIT-TRANS-RXN-L2R', 'SUCROSE-TRANS-RXN-L2R'), n.rxns=c('AMMONIUM-TRANS-RXN-L2R', 'NITRATE-TRANS-RXN-L2R'), s.rxns='SULFATE-TRANS-RXN-L2R')
#rad.rxns <- c(nut.rxns, rad.rxns)
if (fba.lp$obj>10^-6){
    ko.names <- rad0$SYMBOL[match(names(rad.rxns), rownames(rad0))]
    ckr <- call.check.ko(s.ck, rad.rxns, ko.names=ko.names, set.name='Rad', obs=0, annot.df=rad0[names(rad.rxns),], ctrl=list(trace=0, method=1), cprim=numeric(ncol(s.ck)), ub=fva.ub)
}

##nutrients
cat('\n Checking for growth w/o source of: o2, c, n, p, s \n')
nut.rxns <- list(o2='OXYGEN-MOLECULE-TRANS-RXN-L2R', c.rxns=c('CIT-TRANS-RXN-L2R', 'SUCROSE-TRANS-RXN-L2R'), n.rxns=c('AMMONIUM-TRANS-RXN-L2R', 'NITRATE-TRANS-RXN-L2R'), 
p.rxns='Pi-TRANS-RXN-L2R', s.rxns='SULFATE-TRANS-RXN-L2R')
ckn <- check.ko(s.test=s.ck, ko.lst=nut.rxns, sense='E', quiet=FALSE)
cat('checking for growth on D-methionine \n')
cn.dmet <- check.new.media(s.ck, ub=fva.ub, rm='SULFATE', add='CPD-218[CCO-EXTRACELLULAR]', ctrl=list(trace=0), quiet=FALSE)
cat('checking for growth on acetate \n')
cn.a <- check.new.media(s.ck, rm='SUCROSE', add='ACET[CCO-EXTRACELLULAR]', supp.ub=10, ctrl=list(trace=0), quiet=FALSE)

##supplements
cat('\n Testing supplements \n')
rad.supp <- rad0[names(rad.rxns),]
rad.supp <- rad.supp[rad.supp$SUPPLEMENTS!='',]
tes <- test.supp(sp=s.ck, ko.lst=rad.rxns[rownames(rad.supp)], annot=rad.supp, ub=fva.ub, supp.ub=10, cprim=numeric(ncol(s.ck)))

cat('Summary per condition\n')
supp.v <- unlist(tes$supp)
supp.v <- supp.v[!is.na(supp.v)]
cat('\n We get', sum(supp.v>10**-2), 'of n =', length(supp.v), 'conditions correct. \n')
cat('The proportion of supplement conditions we get right', round(mean(supp.v>10**-2), 4), '\n')
cat('The conditions we get wrong are:\n')
print(names(supp.v)[supp.v<10**-2])

cat('Summary per gene\n')
cat('\n Supplement mutants we get wrong \n')
#remove those w/ NA in 2nd row, since these grew w/o supp
tes2 <- tes$mat[,!is.na(tes$mat[2,])]
tes.wrong <- tes2[,tes2[1,]<10**-2 & tes2[2,]<10**-2]
print(rad0[colnames(tes.wrong),c('SYMBOL', 'SUPPLEMENTS.for.Growth', 'ReplaceCPD', 'ReplaceWithCPD')])
print(summary(as.factor(tes2[2,]>10**-2 & tes2[1,]<10**-2)))
