##jmd
##2.24.12
##farm_init.r

source('/msc/neurospora/FBA/farm/farm_config.r')

mm.ss <- get.met.mat(univ.smm=smm.meta, org.smm=smm.nc, s.rownames=rownames(s.sp), chem.form.col='CHEMICAL.FORMULA', colsums.mm.thresh=2)
#write.csv(mm.ss, 'met_mat.csv')
#mm.ss <- as.matrix(read.csv('met_mat.csv', row.names=1))

bal.mat <- get.imbalance(sp=s.sp, met.mat=mm.ss, biomass.rxns=bm$rxn, trans.regexp='-TRANS-RXN')
#need to write this for make_rxn_annot
write.csv(bal.mat, 'rxns_imbalance.csv')

#know from henry (2007) table 4: nad/nadh=15; nadp/nadph=1.1; (atp+.5*adp)/(atp+adp+amp)=0.9 [andersen & meyerburg (jbc, 77)]
#using bennet (nat chem bio, 09) also (in Molar): [nad]=2.5*10^-3 M; [nadh]=[nad]/15; [nadph]=2.5*10^-4; [nadp]=1.1*[nadph]; [amp]=2*10^-4; [adp]=3*10^-4; [atp]=0.003
#water~=1 M, h+[intra]=10**-7 M, o2<10**-5 M, co2 diff than other met concentrations
known.mets.conc=c(PROTON=10**-7, NAD=2.5*10^-3, NADH=(2.5*10^-3)/15, NADPH=2.5*10^-4, NADP=1.1*(2.5*10^-4), AMP=2*10^-4, ADP=3*10^-4, ATP=0.003)

#rr <- make.rxn.annot(sp=s.sp, univ.ec=meta.ec, org.ec=nc.ec, biomass.rxns=bm$rxn, pwy.rxns=nc.pwy$rxn, suffix.regexp='\\[CCO-.+\\]', gpr=gpr, known.mets.conc=known.mets.conc,
#eqn.col='eqn', kegg.col='Kegg', rxn.col='rxn', ecpred.file='/msc/neurospora/FBA/Neurospora/NC10_CALLGENES_FINAL_2.ecpred', rxns.imb.file='rxns_imbalance.csv', write.file='rxn_annot.txt')
rr <- read.delim('rxn_annot.txt')

##filter
#rm rxns that can't carry flux
max.nv <- max.n.flux(s.sp, se='E', allow.all.trans=TRUE)
#rm rxns that add CNOPS elements
sp.f <- ss.mat(s.sp, which(max.nv>10**-6 & rr$unb.add<=0))
rr.f <- rr[colnames(sp.f),]
gpr.f <- gpr[gpr$RXN %in% colnames(sp.f),]

##get non-dilute cpds
#ignore extracell cpds when do md constraints
rs.spf <- rowSums(abs(sp.f))
names(rs.spf) <- rownames(sp.f)
nondilute.cpds <- unique(gsub('\\[CCO-.+', '', names(rs.spf)[rs.spf>=105]))
nondilute.cpds.regexp <- paste(nondilute.cpds, collapse='|')
nondilute.regexp <- '^(nondilute.cpds.regexp)(\\[CCO-|$)|\\[CCO-EXTRACELLULAR\\]$'

##connectivity constraints
cut.met.mat <- cut.mets(sp.f, sense='E', n.top.mets.rm=300)
#cut.rev.mat <- cut.revs(sp.f, n.top.mets.rm=2000)

fva.ub <- rep(10**3, ncol(sp.f)); names(fva.ub) <- colnames(sp.f)
fva.ub[c('CARBON-DIOXIDE-TRANS-RXN-L2R', 'CIT-TRANS-RXN-L2R', 'SUCROSE-TRANS-RXN-L2R')] <- c(0, 0, 5)
#penalize rxns that use met classes
p <- get.probs(rxn.annot=rr.f, prob.cols=c('ef.prob', 'thermo.prob'), gpr.prob=0.8, pwy.prob=0.9, unb.rm.pen=1, h.unb.add.pen=10, unb.add.pen=100)
p[known.rxns] <- 1
