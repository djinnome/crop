##jmd
##6.12.12
##check_stories.r
##checks any stories/statements made in paper but not in training/test set

source('/msc/neurospora/FBA/farm/farm_header.r')

###ace-2,3,4 rescue
s2 <- changeNuts(s.al, add.nuts='L-CITRULLINE')
f <- FBA(s2, ko='PYRUVDEH-RXN-L2R[CCO-MIT]')

s.trace <- s2
#this rxn uses oxoglutarate to produce GLT, so shouldn't count here
s.trace['GLT', 'ORNITHINE-GLU-AMINOTRANSFORASE-RXN-L2R'] <- 0
trace.mets(s.trace, f$xopt, from='L-CITRULLINE', to='2-KETOGLUTARATE', pwy=nc.pwy)

###arg-14 non-rescue by FBA
#to rescue arg-14, FBA would 1st need arg-14 to be lethal. this could be done by adding 'ACETYL-GLU' to biomass.
ko.rxns <- ncu2rxn(ncu=genes[['arg-14']], gpr=gpr)
f <- check.new.media(s.al, ko=ko.rxns, add='ARG', goal='ACETYL-GLU') #no growth - can't rescue

###ad-5 rescue when allow accumulation of AICAR
ko.rxns <- NCUvector2rxns(genes[['ad-5']], gpr=gpr)
se <- named.vec('E', name=rownames(s.al))
se['AICAR'] <- 'G'
f1 <- check.new.media(s.al, ko=ko.rxns, add='ADENINE') #does not grow
f2 <- check.new.media(s.al, ko=ko.rxns, add='ADENINE', se=se) #grows
