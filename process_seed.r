##jmd
##5.27.2011
##process_seed.r
#takes in: ModelSEED-reactions-db, ModelSEED-compounds
#writes out: ModelSEED-reactions, s0_seed, biomass_seed, vogel2, cpd_map, rxn_map?

setwd('/msc/neurospora/FBA/seed')

source('../farm/sparse_mat_fcns.r')

##edit rxn annotations
rxns <- read.csv('/msc/neurospora/FBA/seed/ModelSEED-reactions-db.csv', as.is=TRUE)
#empty therm
rxns$THERMODYNAMIC.FEASIBILTY[rxns$THERM==''] <- '<>'
#for compatibility w/ calc/excel
rxns$THERMODYNAMIC.FEASIBILTY <- gsub('=', '', rxns$THERM)
#not all rxns split by '<=>'
rxns$EQUATION <- gsub('^=> ', ' <=> ', rxns$EQUATION)
#some rxns are a+b<=>b+a, and are " <=> " in EQUATION column
rxns <- rxns[rxns$EQUATION!=" <=> ",]
#rxn 3-methyl-2-oxobutanoate:ferredoxin oxidoreductase (decarboxylating-CoA-2-methylpropanoylating) had equation in EC column
mofo.ind <- which(rxns$NAME=='3-methyl-2-oxobutanoate:ferredoxin oxidoreductase (decarboxylating-CoA-2-methylpropanoylating)')
rxns$EQUATION[mofo.ind] <- gsub('\\|', '', rxns$EC[mofo.ind])
rxns$EC.NUMBER.S.[mofo.ind] <- ''
#'cdp00251' is a misspelling
rxns$EQUATION <- gsub('cdp00251', 'cpd00251', rxns$EQUATION)
#rxn00826, rxn00644, rxn00087, rxn11118 have '2' instead of '(2)'
no.paren.rxns <- c('rxn00826', 'rxn00644', 'rxn00087', 'rxn11118')
no.paren.ind <- which(rxns$DATA %in% no.paren.rxns)
rxns$EQUATION[no.paren.ind] <- gsub('^2 ', '(2) ', gsub(' 2 ', ' (2) ', rxns$EQUATION[no.paren.ind]))
#rxn14377 uses 'n-1' as coeff, HAS 'n-1cpd...'; so let n=2
rxns$EQUATION[rxns$DATA=='rxn14377'] <- gsub('n-1', '', rxns$EQUATION[rxns$DATA=='rxn14377'])
#get KEGG="" sometimes and NA sometimes
rxns$KEGG[rxns$KEGG==''] <- NA
#write new rxns file
write.csv(rxns, '/msc/neurospora/FBA/seed/ModelSEED-reactions.csv', row.names=FALSE)

##get seed s.sp
#added cpd15595 Ala-Ala to 'ModelSEED-compounds-db.csv' manually
rxns <- read.csv('/msc/neurospora/FBA/seed/ModelSEED-reactions.csv', as.is=TRUE)
s.sp <- sp <- spread2sparse(rxns)
write.table(s.sp, 's0_seed_cpds.txt', sep='\t', row.names=FALSE, quote=FALSE)

##replace cpd names
#added cpd15595 to compounds-db.csv
#have to deal w/ compartmentalization suffix
cpds <- read.csv('/msc/neurospora/FBA/seed/ModelSEED-compounds.csv', as.is=TRUE)
cpds.no.compart <- gsub('\\[[[:alpha:]]\\]$', '', sp$comp) 
compound.names <- cpds$PRIMARY.NAME[match(cpds.no.compart, cpds$DATABASE)]
compart.ind <- grep('\\[.+]$', sp$comp)
compart.names <- gsub('^.+\\[|]$', '', sp$comp[compart.ind])
full.comp.names <- compound.names
full.comp.names[compart.ind] <- paste(full.comp.names[compart.ind], '[', compart.names, ']', sep='')
sp$compound <- full.comp.names
#write sp
write.table(sp, 's0_seed.txt', sep='\t', row.names=FALSE, quote=FALSE)

##replace vogel compounds
#write out some matched compounds, then do rest by hand
vogel <- read.delim('../farm_data/vogel_trans_s.txt')
#many vogel$comp didn't match smm$FRAME
vogel2 <- data.frame(vogel, sd.comp=comb$PRIMARY.NAME[match(vogel$compound, comb$FRAME)])
write.table(vogel2, 'vogel2.txt', sep='\t', row.names=FALSE, quote=FALSE)

##replace biomass compounds
#write out some matched compounds, then do rest by hand
bm <- read.delim('../farm_data/biomass_s.txt')
bm.match0 <- data.frame(bm.cpd=unique(bm$comp), cpd.match=comb$PRIMARY.NAME[match(unique(bm$comp), comb$FRAME)], manual=FALSE, notes=NA)
#write.csv(bm.match[order(bm.match$cpd, bm.match$bm),], 'bm_match0.csv', row.names=FALSE, quote=FALSE, na='')
#match by hand then re-read
#swapped 'N-Acetyl-D-glucosamine' for 'chitobiose', since it couldn't get produced organically
#removed beta-Carotene, since can't make it w/o cycles
bm.match <- read.csv('bm_match.csv', as.is=TRUE, row.names=NULL)
bm2 <- bm
bm2$compound <- bm.match$cpd.match[match(bm$comp, bm.match$bm.cp)]
#get rid of unmapped compounds from biomass
bm2 <- bm2[bm2$comp!='',]
#change coeffs
bm2$coeff[bm2$rxn=='biomass' & bm2$comp=='Carbohydrates'] <- -0.1
bm2$coeff[bm2$coeff < -1] <- -1
#write
write.table(bm2, 'biomass_seed.txt', sep='\t', row.names=FALSE, quote=FALSE)

##merge annotations of compounds
#read dbs
cpds <- read.csv('/msc/neurospora/FBA/seed/ModelSEED-compounds.csv', as.is=TRUE)
smm <- read.delim('/msc/neurospora/FBA/Meta/BiomassComposition/metacyc.smm', as.is=TRUE)
#merge, even tho cpds has multi KEGG IDs per line, separated by "|"
smm.int <- smm[smm$KEGG %in% unlist(strsplit(cpds$KEGG, split='|', fixed=TRUE)),]
cpd.int.ind <- unlist(apply(X=smm.int, MARGIN=1, FUN=function(x){ grep(x=cpds$KEGG, pattern=x['KEGG'])[1] }))
comb.cpds <- cbind(cpds[cpd.int.ind, c('DATABASE','PRIMARY.NAME','FORMULA','KEGG.ID.S.')], kegg=smm.int$KEGG, smm.int[,c('FRAME','COMMON.NAME','CHEMICAL.FORMULA')])
#write
write.csv(comb.cpds, 'cpd_map.csv', row.names=FALSE)
