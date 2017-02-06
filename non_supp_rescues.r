##jmd
##2.23.12
##non_supp_rescues.r
##test negative control of mutants x non-supps, but need to validate these non-supps

#source('run_check_ko.r')

map <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/eCompendium/nutrient-biocyc-map.txt')
#get rid of classes, which have some lower case letters
map <- map[-grep('[a-z]', map[,1]),]
map.ss <- map[paste(map[,1], '[CCO-EXTRACELLULAR]', sep='') %in% rownames(s.al),]

rad.supp <- rad0[names(rad.rxns),]
rad.supp <- rad.supp[rad.supp$SUPPLEMENTS!='',]

supps <- parse.supp(annot=rad.supp)

rad.supp2 <- rad.supp
for (i in 1:length(supps)){
    rad.supp2$SUPPLEMENTS.for.Growth[i] <- paste(setdiff(map.ss[,1], unlist(supps[ rownames(rad.supp)[i] ])), collapse=' or ')
}

tt <- test.supp(sp=s.ck, ko.lst=rad.rxns[rownames(rad.supp2)], annot=rad.supp2, sense='E')

cat('Summary per gene\n')
cat('\n Supplement mutants we get wrong \n')
#remove those w/ NA in 2nd row, since these grew w/o supp
tes2 <- tt$mat[,!is.na(tt$mat[2,])]
tes.wrong <- tes2[,tes2[1,]<10**-6 & tes2[2,]>10**-6]
print(rad0[colnames(tes.wrong),c('SYMBOL', 'SUPPLEMENTS.for.Growth', 'ReplaceCPD', 'ReplaceWithCPD')])
print(summary(as.factor(tes2[2,]>10**-6 & tes2[1,]<10**-6)))

cat('Summary per condition\n')
supp.v <- unlist(tt$supp)
supp.v <- supp.v[!is.na(supp.v)]
cat('\n We get', sum(supp.v>10**-6), 'of n =', length(supp.v), 'conditions that grow. \n')
cat('The proportion of supplement conditions for which we get growth', round(mean(supp.v>10**-6), 4), '\n')
cat('The conditions we get growth on are:\n')
supp.wrong <- names(supp.v)[supp.v>10**-6]
print(sort(supp.wrong))
