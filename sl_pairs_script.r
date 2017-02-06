##jmd
##6.5.12
##sl_pairs_script.r

source('/msc/neurospora/FBA/farm/farm_header.r')

non.essential.thresh <- 0.05

#non-essential genes
g.mat <- read.csv('singleKO_growth.csv', row.names=1)
g.mat$growth[is.na(g.mat$growth)] <- 0
ne.genes <- setdiff(genes, g.mat[g.mat$growth<=non.essential.thresh, 'genes'])
ne.genes <- setdiff(ne.genes, genes[c('ace-7', 'ace-8', 'thi-4')])

#construct matrix of gene pairs
sl.gene.mat <- t(combn(x=ne.genes, m=2))
colnames(sl.gene.mat) <- c('g1', 'g2')
rownames(sl.gene.mat) <- apply(sl.gene.mat, 1, FUN=function(x){ paste(x[1], ',', x[2], sep='') })

sl.growth <- multiGeneKO(a=s.al, ncu=sl.gene.mat, gpr=gpr, fba.ub=fva.ub, quiet=TRUE)

length(sl.v <- names(sl.growth)[!is.na(sl.growth) & sl.growth<=10**-7])
names(sl.v) <- sapply(strsplit(x=sl.v, split=','), FUN=function(x){ paste(names(genes)[match(x, genes)], collapse=',') })
sl.v <- sl.v[order(names(sl.v))]

write.csv(sl.v, 'pred_synth_lethals.csv')
