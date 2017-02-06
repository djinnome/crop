##jmd
##8.11.11
##get_sachs_genes.r

setwd('/msc/neurospora/FBA/farm_data')

broad.all <- read.delim('/msc/neurospora/FBA/Neurospora/broad_phenotypes/minimal-media-phenotypes.txt')
rad.all <- read.delim('/msc/neurospora/FBA/Neurospora/eCompendium/genes.txt')
gpr <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.gpr')

done <- union(rad.all$LOCUS, sub('.5', '', broad.all$Locus, fixed=TRUE))
sdf <- setdiff(unique(gpr$Genes), done)
sdf.rxns <- NCUvector2rxns(sdf, gpr)
