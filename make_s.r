##jmd
##2.13.11
##make_s.r

source('/msc/neurospora/FBA/farm/globals.r')
setwd('/msc/neurospora/FBA/farm_data')

##made sparse format biomass_s.txt using make_biomass.r
##created sparse format transport rxns (import Vogel's media) in excel

##read
bm <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.biomass')
bm.goals <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.biomass-goals')
vogel <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.nutrients')
model.rxns <- read.csv('model_rxns.csv')[,1]
nc <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.S')
meta <- read.delim('/msc/neurospora/FBA/farm_data/Meta/meta.S')

##combine s mats
#no dup rxns from meta & nc10
s <- merge.s(univ=meta, org=nc)

##sink
#don't need to send it vogel's, since it already acts on every extracellular metabolite
#exch.df <- make.sink(s=s)
exch.df <- read.delim('/msc/neurospora/FBA/farm_data/sink.tsv')
#setdiff w/ vogels
export <- exch.df[!(exch.df$rxn %in% vogel$rxn),]

#add transport rxns
s <- rbind(s, bm, bm.goals, vogel, export)

##write
s <- s[order(s$rxn),]
write.table(s, 's15.txt', sep='\t', row.names=FALSE, quote=FALSE)

###check#########################################################################################
##check s$rxn
#sort(unique(s$rxn))[1:50]

##unaccounted for transport rxns
#import
imp <- tapply(X=as.numeric(s$coeff), INDEX=s$rxn, FUN=function(x){ sum(x<0)==0 })
ii <- names(imp)[imp]
cat("Number of reactant-less internal rxns:", length(ii <- setdiff(ii,  vogel$rxn)), '\n')
#export
ex <- tapply(X=s$coeff, INDEX=s$rxn, FUN=function(x){ sum(x>0)==0 })
ee <- names(ex)[ex]
cat("Number of product-less internal rxns:", length(ee <- setdiff(ee, c(bm$rxn, bm.goals$rxn, vogel$rxn, export$rxn))), '\n')

##model.rxns
cat('model rxns not in S15 are:', setdiff(model.rxns, s$rxn), '\n')
cat('model rxns in S15 but not in (nc10.s or nc10.export):', setdiff(intersect(model.rxns, s$rxn), c(nc$rxn, bm$rxn, vogel$rxn, export$rxn)), '\n')

#delete from s
#s <- s[!(s$rxn %in% ii),]
#s[s$rxn %in% names(ex)[ex],]
##write
#si <- s[s$rxn %in% ii,]
#write.csv(si, 'free_prod_rxns.csv')
