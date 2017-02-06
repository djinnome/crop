##jmd
##12.6.11
##farm_header.r

source('/msc/neurospora/FBA/farm/farm_config.r')
source('/msc/neurospora/FBA/farm/phenos2rxns.r')

##read
#rxns annot
rxns <- read.delim('rxn_annot.txt')
rxns <- rxns[colnames(s.sp),]
all(rownames(rxns)==colnames(s.sp))

##subset
model.rxns <- read.csv('model_rxns.csv')[,1]
s.al <- ss.mat(sg, c=intersect(model.rxns, colnames(sg)))
rxns.al <- rxns[colnames(s.al),]
nmets <- nrow(s.al); nrxns <- ncol(s.al)

##matched metabolite mass
smm2 <- smm[rownames(s.al),]
met.mass <- named.vec(1, rownames(s.al))
met.mass[is.na(smm2$MASS)] <- 10**5
met.mass[!is.na(smm2$MASS)] <- smm2$MASS[!is.na(smm2$MASS)]

##fba
sp.al <- s.sp[rownames(s.al), colnames(s.al)]
#sp.al[c('NADH', 'FMN'), 'biomass'] <- -0.01

##gene names
gpr <- gpr[gpr$RXN %in% colnames(s.al),]
#all genes
genes <- setdiff(unique(gpr$Genes), 's0001')
g.names <- ncu.gene.map$Symbol[match(gsub('\\.5$', '', genes), ncu.gene.map$Locus)]
g.names[is.na(g.names)] <- genes[is.na(g.names)]
names(genes) <- g.names

##ub
#DO *NOT* change this ub to Inf, it is used in al.a
fva.ub <- rep(10**3, nrxns); names(fva.ub) <- colnames(s.al)
fva.ub['SUCROSE-TRANS-RXN-L2R'] <- 5
#fba
cat('Checking growth:', '\n')
f <- FBA(s.al, sense='E', fba.ub=fva.ub)
f.full <- FBA(s.al, sense='E', fba.ub=fva.ub, fba.obj='FullBiomassComposition')
