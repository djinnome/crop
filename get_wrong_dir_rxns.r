##jmd
##2.13.11
##get_wrong_dir_rxns.r
#finds non-pwy rxns that are reversible in nccyc but not in metacyc

options(stringsAsFactors=FALSE)
source('/msc/neurospora/FBA/farm/globals.r')
setwd('/msc/neurospora/FBA/farm_data')

##made sparse format biomass_s.txt using make_biomass.r
##created sparse format transport rxns (import Vogel's media) in excel

##read
#nc
nc.ec <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.EC')
rownames(nc.ec) <- nc.ec$rxn
nc.ec$dir.rxn <- 'L2R'; nc.ec$dir.rxn[grep('-R2L$', nc.ec$rxn)] <- 'R2L'
nc.pwy <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/WorkingLipids/nc10.pwy')
#meta
meta.ec <- read.delim('/msc/neurospora/FBA/farm_data/Meta/meta.EC')
meta.ec$dir.rxn <- 'L2R'; meta.ec$dir.rxn[grep('-R2L$', meta.ec$rxn)] <- 'R2L'
#mm=meta.ec[meta.ec$rxn %in% meta.ec$rxn[duplicated(meta.ec$rxn)],]; mm <- mm[order(mm$rxn),]
meta.ec <- meta.ec[!duplicated(meta.ec),]
rownames(meta.ec) <- meta.ec$rxn

##get wrong dir rxns in setdiff(nc, nc.pwy)
nc2 <- nc.ec[!(nc.ec$rxn %in% nc.pwy$rxn) & (nc.ec$FRAME %in% meta.ec$FRAME),]
inds <- which(apply(nc2, 1, FUN=function(x){ !(x['dir.rxn'] %in% meta.ec$dir.rxn[meta.ec$FRAME %in% x['FRAME.ID']]) }))
nc.wrong.dir <- nc2[inds,]
#write
nc.wrong.dir <- nc.wrong.dir[order(nc.wrong.dir$rxn),]
write.table(nc.wrong.dir, 'nc_wrong_dir_rxns.txt', sep='\t', quote=FALSE, row.names=FALSE)
