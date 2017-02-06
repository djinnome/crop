##jmd
##8.22.2011
##farm_config.r

##we don't need boudnedness b/c we use goal programming

source('/msc/neurospora/FBA/farm/globals.r')

setwd('/msc/neurospora/FBA/farm_data')

bm <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.biomass')
bm.goals <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.biomass-goals')
vogel <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.nutrients')
nc <- read.delim('Neurospora/nc10.S')
nc.pwy <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.pwy')
gpr <- read.delim('Neurospora/nc10.gpr')
filt.rxns <- read.delim('Neurospora/nc10.filtered'); rownames(filt.rxns) <- filt.rxns$rxn
wrong.dir.rxns <- read.delim('nc_wrong_dir_rxns.txt'); rownames(wrong.dir.rxns) <- wrong.dir.rxns$rxn
model.rxns <- read.csv('/msc/neurospora/FBA/farm_data/model_rxns.csv')[,1]

#gene map has g2t-7 and g2t-8
ncu.gene.map <- read.table('/msc/neurospora/FBA/farm_data/Neurospora/nc10_symbols.txt', sep='\t', header=TRUE, as.is=TRUE)
non.ncu.in.map <- grep('G2T-', ncu.gene.map$Locus)
rownames(ncu.gene.map) <- paste(ncu.gene.map$Locus, '5', sep='.')
rownames(ncu.gene.map)[non.ncu.in.map] <- ncu.gene.map$Locus[non.ncu.in.map]

nc.ec <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.EC')
meta.ec <- read.delim('/msc/neurospora/FBA/farm_data/Meta/meta.EC')

smm.nc <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.smm', na='NIL')
smm.nc <- smm.nc[,!(colnames(smm.nc) %in% c('CPD.ID', 'LOCATION', 'InChI'))]
smm.meta <- read.delim('/msc/neurospora/FBA/farm_data/Meta/meta.smm', na='NIL')

nut.rxns <- list(o2='OXYGEN-MOLECULE-TRANS-RXN-L2R', c.rxns='SUCROSE-TRANS-RXN-L2R', n.rxns=c('AMMONIUM-TRANS-RXN-L2R', 'NITRATE-TRANS-RXN-L2R'), 
p.rxns='Pi-TRANS-RXN-L2R', s.rxns='SULFATE-TRANS-RXN-L2R')
known.rxns <- c(unlist(nut.rxns), bm='biomass', full.bm='FullBiomassComposition')

smm <- read.csv('smm.csv')
rownames(smm) <- smm$FRAME

##s
s.sp <- read.S('s15.txt')
nmets <- nrow(s.sp); nrxns <- ncol(s.sp)

##metab dilute s
#mets.nocco <- sub('\\[CCO-.+$', '', rownames(s.sp))
#snt <- apply(s.sp, 2, FUN=function(x){ tapply(x, INDEX=mets.nocco, FUN=sum) })
#remove pure transport + biomass/full-biomass composition rxns
trans.rxns <- nc.ec$rxn[nc.ec$Transporter.p=='Transporter']
sg <- metab.dilute.s(s.sp, eps=10**-4, nondilute.rxns=c(trans.rxns, unique(bm$rxn)),
compartments.dilute=c('', paste('[CCO-', c('MIT', 'GLYOXYSOME', 'NUC-LUM'), ']', sep='')) )
