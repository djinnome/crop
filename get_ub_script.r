##jmd
##9.13.11
##get_ub_script.r
##code to get tight bounds using fcns

##to do: 
#1. loop st start w/ low default dual.ub's for dual.ub.milp, then set known dual.ub's to new values, and try higher default dual.ub, etc
#2. loop st when get associated dual rxns for one rxn, don't run dual.ub.milp on these but just use FBA w/ v_i=1, w/ rest of dual rxn set KO'd

##run farm_script to get s.al w/o no.v rxns

need.rxns <- read.table('necessary_rxns.txt', header=TRUE)[,1]
#dub <- dual.ub(s.al, need.rxns)
#dub.nona <- dub[!is.na(dub)]
#write.csv(dub.nona, 'dual_ub_nec_rxns.csv')
lub <- rep(10**3, ncol(s.al))
names(lub) <- colnames(s.al)
#need.rxns must be in net, so their beta=1 -> thier lambda <= 0
lub[intersect(need.rxns, names(lub))] <- 0

#md <- 1/dub.nona
#mf <- min.flux(s.al, rxns=need.rxns, min.bm=1)
#all(abs(mf[names(md)]-md)<10**-1) #TRUE
#as fva bounds, these simultaneously do not allow growth of 1.

##dum
f <- FBA(s.al, fba.ub=rep(10**4, ncol(s.al)))
cut.met.mat <- cut.mets(s.al, sense='E')
#test by comparing 1st few to dual.ub
(dumi <- dual.ub.milp(a=s.al, r=need.rxns[1:3]))
(dumi2 <- dum.cut(a=s.al, r=need.rxns[1:3], n.ko=Inf, cut.met.mat=cut.met.mat, ctrl=list(tilim=30)))

#other rxns
noneed <- setdiff(colnames(s.al), need.rxns)
(dumi <- dual.ub.milp(a=s.al, r=noneed[1:3], n.ko=5, lam.ub=lub, fba.ub=rep(10**4, ncol(s.al)), ctrl=list(tilim=30)))
(dumi2 <- dum.cut(a=s.al, r=noneed[1:5], n.ko=3, lam.ub=lub, fba.ub=rep(10**4, ncol(s.al)), ctrl=list(tilim=10), cut.met.mat=cut.met.mat))

(dumi2 <- dum.cut(a=s.al, r='RXN66-312-CPD-4576/NADPH/OXYGEN-MOLECULE//CPD-4577/NADP/WATER.52.-L2R', n.ko=Inf, lam.ub=lub, fba.ub=rep(10**4, ncol(s.al)), ctrl=list(tilim=60), cut.met.mat=cut.met.mat))
(dumi2 <- dual.ub.milp(a=s.al, r='RXN66-312-CPD-4576/NADPH/OXYGEN-MOLECULE//CPD-4577/NADP/WATER.52.-L2R', n.ko=3, lam.ub=lub, fba.ub=rep(10**4, ncol(s.al)), ctrl=list(tilim=10)))

(dumi <- dual.ub.milp(a=s.al, r="PEPDEPHOS-RXN-R2L", n.ko=Inf, fba.ub=rep(10**4, ncol(s.al)), lam.ub=lub, ctrl=list(tilim=600)))
(dumi2 <- dum.cut(a=s.al, r="PEPDEPHOS-RXN-R2L", n.ko=Inf, fba.ub=rep(10**4, ncol(s.al)), lam.ub=lub, ctrl=list(tilim=600), cut.met.mat=cut.met.mat))

#results: {(1.2.1.31-RXN-ALLYSINE/NADP/WATER//CPD-468/NADPH/PROTON.42.-R2L, 21.35), (1.3.1.70-RXN-R2L, 9.58), (RXN66-317-CPD-4580/NADH/OXYGEN-MOLECULE//CPD-4702/NAD/WATER.50.-L2R, 9.73)}

###loop thru rxns
dub.v <- rep(NA, ncol(s.al)); names(dub.v) <- colnames(s.al); dub.v[need.rxns] <- 0
r2g <- rxns2go <- names(dub.v)[is.na(dub.v)]; lr2g <- length(rxns2go)
fba.ub.inf <- rep(Inf, ncol(s.al)); names(fba.ub.inf) <- colnames(s.al)
#highest dual.ub capped at fba.max. if get higher ones from dual.u.nec(), then maybe it's not accounting for vogel's bounds=10**3, in FBA().
f <- FBA(s.al); fba.max <- f$obj
while (length(rxns2go) >= lr2g-25){
    (dumi.tmp <- dual.ub.milp(a=s.al, r=rxns2go[1], n.ko=2, lam.ub=lub, fba.ub=rep(10**4, ncol(s.al)), ctrl=list(tilim=10)))
    #if it works, set this as dual ub, get other dual ubs from lst, & delete from rxns2go
    if (dumi.tmp$mat['stat',1] %in% c(101, 102,107) & abs(dumi.tmp$mat['milp.obj',1]-dumi.tmp$mat['fba.obj',1])<0.01){
        dub.v[rxns2go[1]] <- lub[rxns2go[1]] <- dumi.tmp$mat['fba.obj',1]
        #use fba w/ v_j=1 and ko rest of rxn set to get dual.ub_j; could also use lambda value of rxn_j from dual.ub.milp() if it returned this.
        if (length(dumi.tmp$lst[[1]])>=1){
            for (r2 in dumi.tmp$lst[[1]]){
                fba.ub2 <- fba.ub.inf; fba.ub2[r2] <- 1
                kos <- c(rxns2go[1], setdiff(dumi.tmp$lst[[1]], r2))
                fba.tmp <- FBA(s.al, fba.ub=fba.ub2, ko=kos)
                #if this soln is ok, assign it
                if (fba.tmp$obj<fba.max-10**-5){ dub.v[r2] <- lub[r2] <- fba.tmp$obj }
            }#end for
        }#end if length lst
        #delete current rxn from rxns2go & potentially any others solved in this iteration
    }#end if soln is good
    rxns2go <- setdiff(rxns2go, union(rxns2go[1], names(dub.v)[!is.na(dub.v)]))
}#end while length rxns2go
                
#########################################################################################################################
####try to use lp methods to get dual.ub's for unnecessary rxns, but these are ineffective heuristics to np-hard problem
#########################################################################################################################
##exploit ub's that are fixed due to dual.ub (ie all inputs or all outputs in dual.ub)
#158
fixed.mets <- rownames(s.al)[apply(s.al, 1, FUN=function(x){ all(names(x)[x<0] %in% names(dub)) | all(names(x)[x>0] %in% names(dub)) })]
#98
fixed.rxns <- colnames(s.al)[apply(s.al, 2, FUN=function(x){ all(names(x)[x<0] %in% fixed.mets) | all(names(x)[x>0] %in% fixed.mets) })]
length(setdiff(fixed.rxns, names(dub)))

#max v-eps*v_j st v_bm=1, v[ex]=0
#min v st v_bm=1
#but max may take ineffient route to biomass, & min's are too high
    trans.rxns <- grep('^TRANS-RXN.+(CCO-PM-FUNGI|CCO-EXTRACELLULAR)', colnames(s.al), value=TRUE)
    ssm <- ss.mat(s.al, c=trans.rxns)
    ex.cpds <- setdiff(grep('CCO-EXTRACELLULAR', rownames(ssm), value=TRUE), c('WATER[CCO-EXTRACELLULAR]', 'PROTON[CCO-EXTRACELLULAR]', 'CARBON-DIOXIDE[CCO-EXTRACELLULAR]'))
    ex.rxns <- colnames(ssm)[colSums(ssm[ex.cpds,]>0)>0]
    
    no.rxns <- grep('TRANS-RXN2T-51', colnames(s.al), v=T)
    ff <- FBA(s.al, ko=no.rxns)
    v <- ff$xopt; 
    vex <- v[ex.rxns]; (vex <- vex[vex>0])
        
    fba <- FBA(s.al, ko=ex.rxns)

rr=ret[!is.na(ret)]
