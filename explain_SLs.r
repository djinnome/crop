##jmd
##10.28.12
##explain_SLs.r

#run sl_fromMlab_heat.r until 'get mat' section

source('../farm/explain_SLs_fcns.r')

#ex=wrapMinExport(s.al, c('NCU02482.5','NCU06606.5'), fva.ub)
ex=wrapMinExport(s.al, c('NCU02482.5','NCU04899.5'), fva.ub)
gr.tca15 <- testGeneOfSLwExport(s.al, 'tca-15', c('CPD-85', 'PROPIONYL-COA'), pwys.sl)

cluster.genes <- list(g1=getAssocGenes(pwys.sl, 'tca-15'), g2=getAssocGenes(pwys.sl, 'tca-2'))

gr <- named.vec(NA, cluster.genes$g2)
for (gene in names(gr)){
    gr.tmp <- testGeneOfSLwExport(s.al, gene, 'PYRUVATE[CCO-MIT]', pwys.sl)>0.01
    gr[gene] <- mean(gr.tmp)
}

## classify genes
u <- setdiff(unique(unlist(strsplit(names(cluster.genes[[2]]), split=','))), 'NCU02482.5')
e <- sort(summary(as.factor((gpr[gpr$Genes %in% u, 2]))))

### look at results #####################################################################################################
g1.propionyl <- gr
g1.pyr <- gr; g1.nadh <- gr
mean(g1.nadh) # & mean(g1.pyr) are 94.5%

#coq-2: fadh2, FADH2[cco-glyoxysome], pyr[mit], nadh[mit]
#tca-2: pyr[mit] or nadh[mit] explains all, except suc & tca-15; c('CPD-85', 'PROPIONYL-COA') explains all.
# cluster: pyr[mit] or nadh[mit] explains all, except only 80% of suc, and none of tca-15.

### really explain ######################################################################################################
genes <- c('cys-13', 'cys-14')
(ncu.v <- c(setdiff(genes, ncu.gene.map[,2]), paste(ncu.gene.map[ncu.gene.map[,2] %in% genes, 1], '.5', sep='')))
(ko <- ncu2rxn(ncu.v, gpr=gpr))
s2 <- s.al[,setdiff(colnames(s.al), ko)]
f <- FBA(s2, fba.ub=fva.ub[colnames(s2)])
stopifnot(f$obj<10**-3)
fg <- FBA(s2, fba.ub=fva.ub[colnames(s2)], se='G')
if (fg$obj>10**-3){
    mb <- min.export(s2, w=met.mass)
} else {
    mb <- mma.bin(s2, fva=fva.ub[colnames(s2)], import.wts=met.mass, export.wts=rep(10**3, nrow(s.al)))
}
stopifnot(length(mb)==1)
#allow change suggested by mma.bin, & then see what related rxn is essential
names(mb) <- sub('beta_', '', names(mb))
# get rxns of met
rom <- ss.rom(s2, sub('neg_', '', names(mb)))
se <- named.vec('E', rownames(s.al))
if (length(grep('neg_', names(mb)))>0){
    names(mb) <- sub('neg_', '', names(mb))
    se[names(mb)] <- 'G'
    # model wants to get rid of met, so find which essential rxn(s) produce it 
    rom <- colnames(rom)[ rom[names(mb),]>0 ]
} else {
    se[names(mb)] <- 'L'
    rom <- colnames(rom)[ rom[names(mb),]<0 ]
}
# what related rxn is now essential?
now.ess <- NULL
for (r in rom){
    ff <- FBA(s2, fba.ub=fva.ub[colnames(s2)], se=se, ko=r)
    if (ff$obj<10**-3) now.ess <- c(now.ess, r)
}
ss.mat(sp.al, now.ess)
