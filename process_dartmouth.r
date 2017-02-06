##jmd
##2.2.12
##process_dartmouth.r

##cross dartmouth w/ ecompendium
##remove broad & rad training set

source('/msc/neurospora/FBA/farm/farm_header.r')

###dartmouth#############################################################################################################
cat('processing dartmouth KO set\n')
d <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/Dartmouth/ko-of-model.tsv')
d$Strain.type[d$Strain.type=='heterkaryon'] <- 'heterokaryon'
d <- d[d$Strain!='microconidia',]
d <- d[!duplicated(d),]
rownames(d) <- d$Locus
#any(duplicated(d[,1])) #FALSE
#dr <- d[!(d$Locus %in% names(rad.rxns),]; dr[order(dr$Str),]
d2 <- d[!(d$Locus %in% c(broad, rad0$Locus)),]
gpr.model <- gpr[gpr$RXN %in% colnames(s.al),]
d2 <- d2[rownames(d2) %in% gpr.model$Genes,]

###radford ecompendium###################################################################################################
cat('processing radford ecompendium mutant set\n')
mut <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/ko_mutants_w_comments.tsv', row.names='Locus')
mut2 <- mut[setdiff(rownames(mut), c(broad, rad0$Locus)),]
mut2 <- mut2[mut2$PHENO!='',]

#write.table(mut, 'Neurospora/ko_mutants_f.tsv', sep='\t', row.names=FALSE, quote=FALSE)

###test small essentiality set###########################################################################################
cat('test validation set\n')
val.genes <- list(ad3A='NCU03166.5', ad3B='NCU03194.5', inv='NCU04265.5', leu2='NCU04385.5', pan2='NCU10048.5',
ad2='NCU00177.5', ad5='NCU02629.5', ad9='NCU00843.5', arg4='NCU10468.5', his2='NCU09320.5', his3='NCU03139.5', his4='NCU06974.5', his5='NCU06360.5', his6='NCU07156.5')
#am.enam2='NCU01195.5,NCU01744.5', pro3.ota='NCU01412.5,NCU00194.5',
val.rxns <- NCUvector2rxns(ncu.v=unlist(val.genes), gpr=gpr, gpr.name='nc10.gpr')

#sp
g <- rep(NA, length(val.rxns))
names(g) <- val.genes
for (i in 1:length(val.rxns)){
    fba.tmp <- FBA(s.al, fba.ub=fva.ub, ko=val.rxns[[i]], quiet=FALSE, control=list(trace=0, method=1))
    g[i] <- fba.tmp$obj
}
#round(g, 2)
cat('we got', sum(g>=10**-6), 'essential gene wrong. Specificy on', length(g), 'essential genes is:', mean(g<=10**-6), '\n') #94%

##se
d3 <- d2[d2$S=='homokaryon',]
homo.rxn.lst <- NCUvector2rxns(d3$Locus, gpr=gpr)
fva.ub <- rep(10**3, nrxns); names(fva.ub) <- colnames(s.al); fva.ub['SUCROSE-TRANS-RXN-L2R'] <- 5
d3$growth <- NA
model.genes.noko <- setdiff(rownames(d3), names(homo.rxn.lst))
d3[model.genes.noko, 'growth'] <- 1
d3[names(homo.rxn.lst), 'growth'] <- check.ko(s.test=s.al, ko.lst=homo.rxn.lst, ub=fva.ub, sense='E', obs=1, quiet.fba=TRUE, ctrl=list(trace=0, method=1))$pred

cat(sum(d3$growth>10**-6, na.rm=TRUE), 'of', sum(!is.na(d3$growth)), 'non-essential genes were predicted correctly.\n This gives a sensitivity of:', mean(d3$growth>10**-6, na.rm=TRUE), '\n') #92%

###test supplements######################################################################################################
test.supp.df <- read.delim('test_supplements.tsv')
test.supp.df <- test.supp.df[,setdiff(colnames(test.supp.df), 'COMMENTS')]
test.supp.df <- data.frame(test.supp.df, ReplaceCPD='', ReplaceWithCPD='', SUPPLEMENTS=test.supp.df$add.cpd)
#gsub('rm.cpd', 'ReplaceCPD', gsub('add.cpd', 'ReplaceWithCPD', colnames(test.supp.df)))
rownames(test.supp.df) <- test.supp.df$Locus
#heterokaryon rxns
cat('\n Testing supplements \n')
het.rxns <- NCUvector2rxns(test.supp.df$Locus, gpr=gpr)
ckh <- check.ko(s.al, ko.lst=het.rxns)
#April 16, 2012 -- for supplements, used to select non-growers from full radford x borkovich intersection 
#test.supp.df <- test.supp.df[ckh[ckh$pred<=0, 'name'], ]
test.supp.df <- test.supp.df[intersect(names(val.rxns), ckh[ckh$pred<=0, 'name']), ]
tes <- test.supp(sp=s.al, ub=fva.ub, supp.ub=10, ko.lst=het.rxns, annot=test.supp.df, sense='E', cprim=numeric(ncol(s.al))) 

supp.eps <- 10**-6

cat('Summary per gene\n')
cat('\n Supplement mutants we get wrong \n')
#remove those w/ NA in 2nd row, since these grew w/o supp
tes2 <- tes$mat[,!is.na(tes$mat[2,])]
tes.wrong <- tes2[,tes2[1,]<supp.eps & tes2[2,]<supp.eps]
print(test.supp.df[colnames(tes.wrong),c('symbol', 'SUPPLEMENTS')])
print(summary(as.factor(tes2[2,]>supp.eps & tes2[1,]<supp.eps)))

cat('Summary per condition\n')
supp.v <- unlist(tes$supp)
supp.v <- supp.v[!is.na(supp.v)]
cat('\n We get', sum(supp.v>supp.eps), 'of n =',  length(supp.v), 'conditions correct. \n')
cat('The proportion of supplement conditions we get right', round(mean(supp.v>supp.eps), 4), '\n')
cat('The conditions we get wrong are:\n')
print(names(supp.v)[supp.v<supp.eps])
