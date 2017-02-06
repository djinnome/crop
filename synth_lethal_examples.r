##jmd
##3.4.12
##synth_lethal_examples.r
#FBA also gets am;en(am)-2 correct

source('/msc/neurospora/FBA/farm/farm_header.r')

sl <- c(pro3.ota='NCU01412.5,NCU00194.5', pro3='NCU01412.5', ota='NCU00194.5', 
pyr1.uc5='NCU06532.5,NCU07334.5', pyr1='NCU06532.5', uc5='NCU07334.5', 
pro3.aga='NCU02333.5,NCU01412.5', 
arg5.12='NCU01667.5,NCU05410.5', arg5='NCU01667.5', arg12='NCU05410.5',
am.enam2='NCU01195.5,NCU01744.5', am='NCU01195.5', enam2='NCU01744.5')

#call.check.ko() automatically adds [CCO-EXTRACELLULAR]
annot.df <- data.frame(SYMBOL=names(sl), ReplaceCPD=rep('', length(sl)), 
ReplaceWithCPD=c('PRO', '', '', '', rep('URACIL', 3), rep('L-ORNITHINE', 3), rep('', 3)))
rownames(annot.df) <- unlist(sl)
ko.rxns <- NCUvector2rxns(ncu.v=sl, gpr=gpr, gpr.name='model.gpr')

##with supplements
ckk <- call.check.ko(sp=s.al, ko.rxns=ko.rxns, set.name='synth.lethals', sense='E', obs=0, annot.df=annot.df)
rownames(ckk) <- names(sl); ckk

##w/o supps
annot.df$ReplaceWithCPD <- ''
ckk2 <- call.check.ko(sp=s.al, ko.rxns=ko.rxns, set.name='synth.lethals', sense='E', obs=0, annot.df=annot.df)
rownames(ckk2) <- names(sl); ckk2
