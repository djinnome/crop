##jmd
##3.5.12
##run_biolog_confusion.r

source('/msc/neurospora/FBA/farm/globals.r')

pb <- run.pred.biolog(add.trans=TRUE, plate.names='PM1', exclude.trans.mets=c('CIT'), na.is.nogrowth=TRUE, 
extracellular=TRUE, write.out.tsv='biolog_predict', biol.nuts.file='biolog.nutrients', 
biol.wells.file='biolog_wells.tsv', ko.rxns=NULL)

biolog.confusion(growth.low=2, nogrowth.high='1', na.rm=FALSE,
exp.biolog.file=c('pm1-vogels-0.1-day4.tsv', 'pm1-water-0.01-day3.tsv'), 
pred.biolog.file='biolog_predict.tsv', 
levels.ordered=c(0, '<1', 1, '1+', 2, '2+', 3, '3+'))
