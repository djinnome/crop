##jmd
##4.4.11
##get_nut_uptakes.r
##want low ub on nutrient uptakes, to get low fva ub on fluxes

setwd('/msc/neurospora/FBA/farm')

source('farm_header.r')

##use asilomar model: later we'll use experimental data on metacyc 15.0
ra <- read.delim('rxns_asilomar.txt', as.is=TRUE)
rownames(ra) <- ra$rxn
#sa need not match ra, since only use sa from here on
sa <- s.sp[,unique(c(ra$rxn, vogel$rxn))]

##run fba to get max growth
fba.obj <- fba.lb <- numeric(ncol(sa))
fba.ub <- rep(1000,ncol(sa))
names(fba.obj) <- names(fba.ub) <- names(fba.lb) <- colnames(sa)
b <- numeric(nrow(sa))
fba.ub[c('SUCROSE-TRANS-RXN-L2R','CIT-TRANS-RXN-L2R')] <- 20
fba.obj['biomass'] <- 1
#lp
fba.lp <- Rcplex(cvec=fba.obj, Amat=sa, bvec=b, lb=fba.lb, ub=fba.ub, sense='G', objsense='max')
fba.lp$obj #1

##min nut uptakes st high growth
min.nut.obj <- min.nut.lb <- fba.lb
min.nut.lb['biomass'] <- 1
#see if i can set all nut ups<=20
min.nut.lb[vogel$rxn] <- 20
min.nut.obj[vogel$rxn] <- 1
#lp
mn.lp <- Rcplex(cvec=min.nut.obj, Amat=sa, bvec=b, lb=min.nut.lb, ub=fba.ub, sense='G')
y <- mn.lp$xopt; names(y) <- names(fba.ub)

##verify that these ub's give fba growth
ub2 <- fba.ub; mn.v <- ub2[vogel$rxn] <- y[vogel$rxn]
#lp
fba2 <- Rcplex(cvec=fba.obj, Amat=sa, bvec=b, lb=fba.lb, ub=ub2, sense='G', objsense='max')
fba2$stat #1
fba2$obj #1

##write
#need 455 ammonium & 130 Pi, can set everything else to 100
write.csv(mn.v, 'min_nut_uptake.csv')

