##jmd
##4.13.12
##all_supps_by_genes.r

organize.gene.supp.mat <- function(gsm, smm){    
    #gsm[gsm>=4] <- gsm[gsm>=4]-3
    
    #hc <- hclust(dist(t(gsm2)))
    col.o <- c('SULFATE', 'SO3', 'S2O3', 'HOMO-CYS', 'L-CYSTATHIONINE', 'CYS', 'MET', 'MET, THR', 
    'THR', 'ARG', 'L-CITRULLINE', 'L-ORNITHINE', 'LYS', 'LEU', 'HIS', 'INDOLE', 'TRP', 'PHE, TYR', 'TYR', 'ASN', 'HOMO-SER', 'GLN', 
    'ADENINE', 'HYPOXANTHINE',      'NITRITE', 'AMMONIA', '3-HYDROXY-L-KYNURENINE',     'ACET', 'SUCROSE', 'GLC', 'BETA-D-FRUCTOSE', 
    'PANTOTHENATE', 'RIBOFLAVIN', 'NIACINE', 'NIACINAMIDE', '3-HYDROXY-ANTHRANILATE')
    
    gene.supp.mat.o <- gsm[,intersect(col.o, colnames(gsm))]
    #fix metabolite labels
    colnames(gene.supp.mat.o)[-grep(', ', colnames(gene.supp.mat.o))] <- smm[colnames(gene.supp.mat.o)[-grep(', ', colnames(gene.supp.mat.o))], 'COMMON.NAME']
    colnames(gene.supp.mat.o)[colnames(gene.supp.mat.o)=='MET, THR'] <- paste(smm[c('MET', 'THR'), 'COMMON.NAME'], collapse='+')
    colnames(gene.supp.mat.o)[colnames(gene.supp.mat.o)=='PHE, TYR'] <- paste(smm[c('PHE', 'TYR'), 'COMMON.NAME'], collapse='+')
    return(gene.supp.mat.o)
}

source('/msc/neurospora/FBA/farm/farm_header.r')

##training set
rad.supp <- rad0[names(rad.rxns),]
rad.supp <- rad.supp[rad.supp$SUPPLEMENTS!='',]
#remove mistakes on essentials
ko.names <- rad0$SYMBOL[match(names(rad.rxns), rownames(rad0))]
ckr <- call.check.ko(s.al, rad.rxns, ko.names=ko.names, set.name='Rad', obs=0, annot.df=rad0[names(rad.rxns),], ctrl=list(trace=0, method=1), ub=fva.ub)
rad.supp <- rad.supp[intersect(ckr$name[ckr$pred<=10**-6], rownames(rad.supp)),]

##test set
test.supp.df <- read.delim('test_supplements.tsv')
test.supp.df <- test.supp.df[,setdiff(colnames(test.supp.df), 'COMMENTS')]
test.supp.df <- data.frame(test.supp.df, ReplaceCPD='', ReplaceWithCPD='', SUPPLEMENTS=test.supp.df$add.cpd)
rownames(test.supp.df) <- test.supp.df$Locus
#essentials
val.genes <- list(ad3A='NCU03166.5', ad3B='NCU03194.5', inv='NCU04265.5', leu2='NCU04385.5', pan2='NCU10048.5', ad2='NCU00177.5', ad5='NCU02629.5', ad9='NCU00843.5', 
arg4='NCU10468.5', his2='NCU09320.5', his3='NCU03139.5', his4='NCU06974.5', his5='NCU06360.5', his6='NCU07156.5')
val.rxns <- NCUvector2rxns(ncu.v=unlist(val.genes), gpr=gpr, gpr.name='nc10.gpr')
#remove mistakes on essentials
ckh <- check.ko(s.al, ko.lst=val.rxns, ctrl=list(trace=0, method=1), ub=fva.ub)
test.supp.df <- test.supp.df[ckh[ckh$pred<=0, 'name'], ]

##combine sets
col.names <- c('Locus', 'symbol', 'ReplaceWithCPD', 'ReplaceCPD', 'SUPPLEMENTS')
rs2 <- rad.supp[, c('Locus', 'SYMBOL', 'ReplaceWithCPD', 'ReplaceCPD', 'SUPPLEMENTS.for.Growth')]
colnames(rs2) <- col.names
tsd <- rbind(rs2, test.supp.df[,col.names])
all.rxns <- c(rad.rxns, val.rxns)
#predict - used to filter NA preds
tes <- test.supp(sp=s.al, ko.lst=all.rxns[rownames(tsd)], annot=tsd, sense='E', ub=fva.ub, supp.ub=10)
#all non-NA preds from test.supp.df are counted as observed
obs.supps <- unlist(tes[[2]])

##make matrix for heatmap
smm <- read.csv('smm.csv'); rownames(smm) <- smm$FRAME
gene.supp.mat <- get.gene.supp.mat(test.supp.df=tsd, sp=s.al, rxns.lst=all.rxns[rownames(tsd)], obs.supp.names=names(obs.supps), 
supp.col='SUPPLEMENTS', rm.cpds=c('CPD2T-61'), ub=fva.ub, supp.ub=10, eps=10**-2)
rownames(gene.supp.mat) <- tsd$symbol
gsm.o <- organize.gene.supp.mat(gene.supp.mat, smm=smm)
write.table(gsm.o, 'docs/Figures/EssentialsAndSupplements/train_test_supps_x_genes.tsv', quote=FALSE, sep='\t')

##plot
#sg <- as.matrix(read.delim('docs/Figures/train_test_supps_x_genes.tsv', row.names=1, as.is=TRUE))
gs <- gsm.o
#split train + test for plot
rows <- list(train=intersect(rownames(gs), rad0$SYMBOL), test=setdiff(rownames(gs), rad0$SYMBOL))
gs2 <- rbind(gs[rows[[1]],], empty=0, gs[rows[[2]],])
empty.row <- which(rownames(gs2)=='empty')
rownames(gs2)[empty.row] <- ''

##plot transpose
color.v <- c('white', 'red', 'blue', 'seagreen')
pdf(file='docs/Figures/EssentialsAndSupplements/train_test_supps_x_genes_transpose.pdf', height=0.75*7)
dens <- 50
z0 <- z <- gs2
z0[z0>=4] <- 0
make.grid(z0, mar=c(9,9,4,2), color.v=color.v, cex.axis=0.4)
#shading lines
my.rectangle(z, equ=4, col=color.v[2], dens=dens)
my.rectangle(z, equ=5, col=color.v[3], dens=dens)
my.rectangle(z, equ=6, col=color.v[4], dens=dens)
#clear lines in "white space" between test + train rows & add borders to it
segments(y0=c(0.5, 1:ncol(gs2)+0.5), x0=empty.row-0.5, x1=empty.row+0.5, col='white')
abline(v=c(empty.row+0.5, empty.row-0.5))
dev.off()

#########################################################################################################################
#### check stories
#########################################################################################################################
gsm <- get.gene.supp.mat(test.supp.df=tsd, sp=s.al, rxns.lst=all.rxns[rownames(tsd)], supp.col='SUPPLEMENTS', obs.supps=obs.supps, rm.cpds='CPD2T-61')
rownames(gsm) <- tsd$symbol
smm.supp <- smm[colnames(gsm),]
## nit-3
smm.supp[grep('N', smm.supp$CHEM),]

##oxD
smm.supp[grep('S', smm.supp$CHEM),]
