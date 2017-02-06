##jmd
##5.24.11
##fva.r
##assumes all v>=0

##PROBLEM: elementwise UBs do not account for cofactors which are recycled and do not need to be produced from nutrients, eg CoA

library('Matrix')
library('Rcplex')
options(stringsAsFactors=FALSE)

setwd('/msc/neurospora/FBA/farm_data')

##read
rxns <- read.delim('rxn_annot.txt')
bm <- read.delim('biomass_s.txt')
vogel <- read.delim('vogel_trans_s.txt')
s0 <- read.delim('s15.txt', as.is=FALSE, na='NIL')
s.sp <- sparseMatrix(i=as.numeric(s0$compound), j=as.numeric(s0$rxn), x=as.real(s0$coeff))
dimnames(s.sp) <- list(levels(s0$compound), levels(s0$rxn))
nmets <- nrow(s.sp); nrxns <- ncol(s.sp)
#code uses s as s.sp
s <- s.sp

##set up
fba.obj <- lb <- numeric(ncol(s))
#make this max(vgl.flow), from below
fva.ub <- rep(360,ncol(s))
names(fba.obj) <- names(fva.ub) <- names(lb) <- colnames(s)
b <- numeric(nrow(s))
#nut uptake ub's from 'get_nut_ups.r'
fva.ub[vogel$rxn] <- 20
#set these ub's slightly higher than min from 'get_nut_ups.r'
fva.ub['AMMONIUM-TRANS-RXN-L2R'] <- 100; fva.ub['Pi-TRANS-RXN-L2R'] <- 30
fba.obj['biomass'] <- 1
#10^-6 was too small for rxn 582, but don't want to KO rxns w/ only small flux. set dynamically?
eps <- 10^-6 

##create elementwise ub's
#get met mat, to account for flow of elements thru rxns
mm <- as.matrix(read.csv('met_mat.csv', row.names=1))
mm <- mm[,c('C','N','O','P','S')]
#separate out different rxn directions
s.pos <- s.neg <- s.sp
s.pos[s.sp<0] <- 0
s.neg[s.sp>0] <- 0
#vogel
#biotin is a cofactor, thioredoxin a protein, so (i think) these are not catabolized
vogel.cat.rxns <- setdiff(vogel$rxn, c("THIOREDOXIN-TRANS-RXN-L2R", "BIOTIN-TRANS-RXN-L2R"))
vgl.met.flow <- t(s.sp[,vogel.cat.rxns]) %*% mm
vgl.flow <- (t(vgl.met.flow) %*% fva.ub[vogel.cat.rxns])[,1]
names(vgl.flow) <- colnames(vgl.met.flow)
#sneaky: use relationship between max & abs to get max at each cell over 2 mats quickly
#rxns.met.flow <- (t(s.pos)%*%mm + t(s.neg)%*%mm + abs(t(s.pos)%*%mm-t(s.neg)%*%mm))/2
#flow is limited by reactants
rxns.met.flow <- t(-s.neg) %*% mm
#new ub
el.ub <- apply(rxns.met.flow, MARGIN=1, FUN=function(x){ min(c(vgl.flow[x>0]/x[x>0], 1000)) })

##combine ub's
comb.ub <- apply(cbind(fva.ub,el.ub), 1, min, na.rm=TRUE)
#write
write.csv(comb.ub, 's_comb_ub.csv')

#performs sv>=eps*|s|*v so don't get v=1000 when there are unfed island cycles
#checks for loss of growth, in case growth relies on unfed cycles which generate new mass
#can permute order to increasing el.ub 
fba.growth <- 0.2
#fva <- function(s, lb, ub, eps=10^(-6)){
    #fva.ind <- 1:ncol(s)
    fva.ind <- which(rxns$nc.pwy==1)
    for (i in fva.ind){
        fva.ub.tmp <- comb.ub
        #find reverse (opposite) of rxn, if exists
        if (rxns$rev[i]==1){ 
            opp.rxn.ind <- setdiff(which(gsub('-L2R$|-R2L$', '', rxns$rxn) %in% gsub('-L2R$|-R2L$', '', rxns$rxn[i])), i)
            fva.ub.tmp[opp.rxn.ind] <- 0
        }
        obj.tmp <- numeric(ncol(s))
        obj.tmp[i] <- 1
        #rhs <- eps*(s[,i]<0)
        ft <- Rcplex(cvec=obj.tmp, Amat=s, bvec=b, ub=fva.ub.tmp, lb=lb, sense='G', objsense='max', control=list(trace=0))
        stopifnot(ft$stat %in% c(1,3))
        #ft$obj
        comb.ub[i] <- ifelse(ft$stat==1 & ft$obj>0, ft$obj, 0)
        #verify effects on fba
        #fba.lp <- Rcplex(cvec=fba.obj, Amat=s, bvec=b, ub=comb.ub, lb=lb, sense='G', objsense='max', control=list(trace=0))
        #stopifnot(fba$obj>0)
        #don't want changes that strongly diminish growth, if so undo (and mark) by making 1001
        #if (fba.lp$obj>0.9*fba.growth){ fba.growth <- fba.lp$obj } else { comb.ub[i] <- 1001 }
        if (i %% 500 == 0) print(paste("Iteration", i))
    }
    jdc(s,c=i)

#test
fba <- Rcplex(cvec=fba.obj, Amat=s.sp, bvec=b, ub=comb.ub, lb=lb, sense='G', objsense='max')
fba$stat
fba$obj

##nc.pwy
ub.pwy <- comb.ub[rxns$nc.pwy==1]
write.csv(ub.pwy[ub.pwy==0], 'no_flux_pwy_rxns.csv')

##write
#still get lots of ub=1000
write.csv(fva.ub, 's_fva_ub.csv')

###validate######################################################################################
#test fba growth
fba <- Rcplex(cvec=fba.obj, Amat=s, bvec=b, ub=comb.ub, lb=lb, sense='G', objsense='max')
fba$stat
fba$obj
x <- fba$xopt; names(x) <- names(lb); x[unique(bm$rxn)]
