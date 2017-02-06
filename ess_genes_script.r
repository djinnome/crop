##jmd
##1.23.12
##ko_all.r

source('/msc/neurospora/FBA/farm/farm_header.r')

essential.genes.thresh <- 10**-6

#each gene's ko'd rxns, which are in model
ko.rxns <- NCUvector2rxns(ncu.v=genes, gpr=gpr, gpr.name='model GPR')
#intitialize matrix of genes that individually KO >= 1 rxn
g.mat <- data.frame(genes=names(ko.rxns), growth=NA)
rownames(g.mat) <- names(genes)[match(g.mat$genes, genes)]
#loop
g.mat$growth <- sapply(ko.rxns, FUN=function(x){
    fba.tmp <- fba.na(a=s.al, ko=x, fba.ub=fva.ub, quiet=TRUE, control=list(trace=0))
    ret <- signif(fba.tmp$obj, 3)
})

ess.ind <- which(g.mat$growth<=essential.genes.thresh)
nrow(essential.genes <- g.mat[ess.ind,])
(essential.genes <- essential.genes[order(rownames(essential.genes)),])

write.csv(essential.genes, 'pred_essentials.csv', quote=FALSE)
write.csv(g.mat, 'singleKO_growth.csv', quote=FALSE)
