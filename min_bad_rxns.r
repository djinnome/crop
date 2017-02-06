##jmd
##6.23.11
##min_bad_rxns.r

write.files <- FALSE
thresh <- 10**-6

#refine this by using eficaz + thermo probs to choose penalties for balanced rxns

source('/msc/neurospora/FBA/farm/farm_config.r')

##read
s.all <- sg
rxns.all <- read.delim('rxn_annot.txt')
rxns.all <- rxns.all[colnames(s.all),]
pen.v <- numeric(nrow(rxns.all)); names(pen.v) <- rownames(rxns.all)
cat('Check that annotations match matrix:', all(colnames(s.all)==rownames(rxns.all)), '\n')
#ub
fva.ub <- rep(10**3, ncol(s.all)); names(fva.ub) <- colnames(s.all)
fva.ub[grep('-TRANS-RXN-L2R$', vogel$rxn, val=TRUE)] <- 1000
trans.rxns <- grep('^TRANS-RXN', colnames(s.all), value=TRUE); fva.ub[trans.rxns] <- 10**6

##penalties
farm.rid <- read.csv('farm_rid_rxns.csv')[,1]
#filt, except prot-mod-rxns, trans-rxns, and CANNOT-BALANCE-RXNS since discriminate against unb.rxns which generate mass
filt.rxns <- filt.rxns[!(filt.rxns$Reason %in% c('PROTEIN-MODIFICATION-RXNS', 'TRANSPORT-RXNS', 'CANNOT-BALANCE-RXNS')),]
pen.v[rownames(rxns.all) %in% filt.rxns$rxn] <- 0.1
#thermo
dg0 <- rxns.all$dg0; names(dg0) <- rownames(rxns.all)
dg0[is.na(dg0)] <- 6; dg0[c(bm$rxn, vogel$rxn)] <- 0
tp <- rxns.all$thermo.prob; names(tp) <- rownames(rxns.all)
tp[is.na(tp)] <- 0.5
#tp <- rxns.all$thermo.prob; tp[is.na(tp)] <- 0; names(tp) <- rownames(rxns.all)
pen.v[rownames(rxns.all) %in% wrong.dir.rxns$rxn] <- 0.1
pen.v[dg0>=20 & rxns.all$nc.pwy==0] <- dg0[dg0>=20 & rxns.all$nc.pwy==0]/10
pen.v[rxns.all$nc.pwy==1 & tp<0.5] <- dg0[rxns.all$nc.pwy==1 & tp<0.5]/10
#meta
pen.v[rxns.all$nc==0] <- pen.v[rxns.all$nc==0]+0.1+5*(1-rxns.all$prob[rxns.all$nc==0])
#farm
pen.v[intersect(farm.rid, names(pen.v))] <- 10
#unbalanced rxns
pen.v[rxns.all$prob<0 & rxns.all$nc==1] <- 10
pen.v[rxns.all$prob<0 & rxns.all$nc==0] <- 20
#naughty rxn
pen.v[intersect(farm.rid, names(pen.v))] <- 10
#3742: uses folylpolyglutamate(n+1); 
naughty.rxns <- c('RXN-3742-R2L', 'NAD+-ADP-RIBOSYLTRANSFERASE-RXN-R2L', 'UPPSYN-RXN-R2L')
pen.v[intersect(names(pen.v), naughty.rxns)] <- 10**3

##if model.rxns in workspace, apply it to pen.v
if ('model.rxns' %in% ls()){
    pen.v[intersect(model.rxns, names(pen.v))] <- 0
    pen.v[setdiff(names(pen.v)[pen.v==0], model.rxns)] <- 0.1
}

##run min.pen
(mp <- min.pen(sp=s.all, goals=c('biomass', 'ISOCIT-CLEAV-RXN-L2R[CCO-GLYOXYSOME]'), pen.v=pen.v))
rxns[rownames(mp),]
(mpm <- min.pen.milp(sp=s.all, pen.v=pen.v, rxns0=rownames(mp), goals=c('biomass', 'ISOCIT-CLEAV-RXN-L2R[CCO-GLYOXYSOME]'), n.add=4, ctrl=list(tilim=180)))
rxns[rownames(mpm),]

##get min adjust
rr <- rxns.all[rownames(mpm),]
rxns1 <- union(rownames(rr), names(pen.v)[pen.v==0])

cat('Number of needed meta rxns:', nrow(meta.add.rxns <- rr[rr$nc==0,]), '\n')
cat('Number of needed filtered rxns:', nrow(filt.add.rxns <- rr[rownames(rr) %in% filt.rxns$rxn,]), '\n')
cat('Number of needed wrong direction rxns:', nrow(wrong.dir.add.rxns <- rr[rownames(rr) %in% wrong.dir.rxns$rxn,]), '\n')
cat('Unbalanced rxns to be added:', unb.add <- intersect(rownames(rxns.all)[rxns.all$prob<0], rxns1), '\n')
cat('Number of bad Gibbs rxns to be added:', length(gibbs.add <- intersect(rownames(rxns.all)[rxns.all$dg0>20], rxns1)), '\n')
cat('farm rid rxns to be added:', farm.add <- intersect(farm.rid, rxns1), '\n')

##check ko
source('/msc/neurospora/FBA/farm/phenos2rxns.r')
#rad
ck <- check.ko(s.test=s.all[,rxns1], ko.lst=rad.rxns, sense='G'); table(ck$obs, ck$pred>10^-6)
#broad
#ckb <- check.ko(s.test=s.all[,model.rxns], rxns.test=rxns.all[model.rxns,], ko.lst=broad.rxns, obs=1, ub=fva.ub[model.rxns]); table(ckb$obs, ckb$pred>10^-6)

##ex
me <- min.export(s.all[,rxns1], v.min=10)

##write
model.rxns <- rxns1
#write.csv(model.rxns, 'model_rxns.csv', row.names=FALSE)
