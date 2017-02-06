##jmd
##6.26.12
##move2matlab.r

###write out supp rescue data############################################################################################
writeSuppData2file <- function(annot){
    supps0 <- unlist(sapply(parse.supp(annot), FUN=function(x){
        sapply(x, FUN=function(y){ 
            if (length(y)>1) y <- paste(y, collapse=' and ')
            return(y) 
        })
    }))
    colnames(annot) <- tolower(colnames(annot))
    supps <- data.frame(locus=gsub('\\.5.', '.5', names(supps0)), supp=supps0)
    supp.df <- data.frame(annot[supps$locus, c('symbol', 'locus', "replacecpd", "replacewithcpd")], supp=supps$supp)
    supp.df[supp.df==''] <- 'NaN'
    write.table(supp.df, '', row.names=FALSE, quote=FALSE, sep='\t')
    return(TRUE)
}

#writeSuppData2file(rad.supp)
#writeSuppData2file(test.supp.df)

# get 'all_supps_by_genes.r'
gsm.o <- read.table('../farm_data/docs/Figures/EssentialsAndSupplements/train_test_supps_x_genes.tsv', sep='\t')
colnames(gsm.o) <- gsub('.', '-', fixed=TRUE, gsub('^X','', colnames(gsm.o)))
colnames(gsm.o)[colnames(gsm.o) %in% smm$COMMON.NAME] <- smm$FRAME[match(colnames(gsm.o)[colnames(gsm.o) %in% smm$COMMON.NAME], smm$COMMON.NAME)]
#can't use '%in%' b/c then lose dimensions
w <- which(gsm.o==2|gsm.o==3|gsm.o==5|gsm.o==6, arr.ind=TRUE)
dd <- data.frame(rownames(gsm.o)[w[,1]], colnames(gsm.o)[w[,2]])
dd <- dd[order(dd[,1], dd[,2]),]
write.table(dd, '', quote=FALSE, row.names=FALSE)

###write out test non-essentials#########################################################################################
gpr <- gpr[gpr$RXN %in% model.rxns,]
test.genes <- intersect(broad, gpr$Gene)

###compare pred ess: r vs mlab###########################################################################################
setwd('C:/Documents and Settings/jdreyf/Desktop/farm_data')
r <- read.csv('singleKO_growth.csv', row.names=1, as.is=TRUE)
r2 <- r[r$growth<0.1,]
mlab <- read.table('matlab_data/allGenesGrowth.txt', as.is=TRUE, header=TRUE)
rownames(mlab) <- mlab$gene
mlab2 <- mlab[is.na(mlab$growth) | mlab$growth<0.1,]

setdiff(r2$genes, mlab2$gene) #nunca! on 27sept12
(setdiff(mlab2$gene, r2$genes))

##check actual ess lists
r.ess <- read.csv('pred_essentials.csv', row.names=1, as.is=TRUE)
mlab.ess <- read.table('matlab_data/allEssGenes.txt', as.is=TRUE)[,1]

setdiff(r.ess$gene, mlab.ess)
mm <- setdiff(mlab.ess, r.ess$gene)
r[r$genes %in% mm,]

#so, new ess genes in mlab all had growth=10e-5 before
#"NCU00712.5" "NCU01175.5" "NCU02305.5" "NCU03633.5" "NCU05982.5" "NCU07719.5" "NCU08671.5" "NCU11381.5" (3 of these are named: coq-1, coq-3, & fpp)

###compare synth lethals#################################################################################################
#continued from above
ml.sl <- read.table('matlab_data/synLethals.txt', as.is=TRUE, header=TRUE)
# look at growth of non-essentials in R
sort(r$growth[r$gene %in% unique(c(ml.sl[,1], ml.sl[,2]))]) #all > 0.3

##on tin
source('/msc/neurospora/FBA/farm/farm_header.r')
ml.sl <- read.table('matlab_data/synLethals.txt', as.is=TRUE, header=TRUE)
sl.growth <- multiGeneKO(a=s.al, ncu=ml.sl, gpr=gpr, fba.ub=fva.ub, quiet=TRUE)
sort(sl.growth) #all < 1.5e-4, most at 10**-5!
