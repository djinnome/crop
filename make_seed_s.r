##jmd
##5.27.11
##make_seed_s.r
#assumes sv>=0

options(stringsAsFactors=FALSE)

setwd('/msc/neurospora/FBA/seed')

##add rxns
bm <- read.delim('biomass_seed.txt')
vogel <- read.delim('vogel_seed.txt')
s0 <- read.delim('s0_seed.txt')
s <- rbind(s0, bm, vogel)

##order
s <- s[order(s$rxn),]

#write
write.table(s, 's_seed.txt', sep='\t', row.names=FALSE, quote=FALSE)
