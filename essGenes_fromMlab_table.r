##jmd
##12.11.12
##essGenes_fromMlab_table.r

source('farm_header.r')

ess <- read.delim('matlab_data/allEssGenes.txt')[,1]

ess.tab <- cbind(Locus=ess, Name=ess)
rownames(ess.tab) <- ess

ess.tab[intersect(ess, rownames(ncu.gene.map)), 'Name'] <- ncu.gene.map[intersect(ess, rownames(ncu.gene.map)), 'Symbols']

ess.tab.o <- ess.tab[order(ess.tab[,'Name']),]

write.table(ess.tab.o, 'matlab_data/allEssGenes_w_names.tsv', row.names=FALSE, quote=FALSE)
