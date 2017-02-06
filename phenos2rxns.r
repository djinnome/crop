##jmd
##6.16.11
##phenos2rxns.r

##ideally this'd be a fcn that takes in model.rxns, gpr, and filenames. maybe standardize input files?

options(stringsAsFactors=FALSE)

gpr <- gpr[gpr$RXN %in% model.rxns,]

##broad
print('Broad KOs')
broad0 <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/broad_phenotypes/minimal-media-phenotypes.txt')
broad <- broad0[,c(1,5)]
broad <- broad[!duplicated(broad),]
#get rid of KOs w/ mult viabilities, leaves only viable phenos
dup.loci <- broad$Locus[duplicated(broad$Locus)]
broad <- broad[!(broad$Locus %in% dup.loci), 'Locus']
#get rxns
print('The following genes are filtered:')
broad.rxns <- NCUvector2rxns(ncu.v=broad, gpr=gpr, gpr.name='model.gpr')

##radford
rad0 <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/eCompendium/essential_gene_phenotypes.txt')
rad0 <- rad0[,colSums(!is.na(rad0))>0]
#name rad rows, but broad NCUs aren't unique
rownames(rad0) <- rad0$Locus
#make rxn list
print('Radford KOs')
print('The following genes are filtered:')
#rad0 <- rad0[order(rad0$EC),]
rad.rxns <- NCUvector2rxns(ncu.v=rad0$Locus, gpr=gpr, gpr.name='model.gpr')
