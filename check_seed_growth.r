##jmd
##5.23.2011
##check_seed_growth.r

library('Matrix')
library('Rcplex')
options(stringsAsFactors=FALSE)

setwd('/msc/neurospora/FBA/seed')

source('../farm/fba.r')

##read
bm <- read.delim('biomass_seed.txt')
vogel <- read.delim('vogel_seed.txt')
s0 <- read.delim('s_seed.txt', as.is=FALSE)
#problems w/ s
s0 <- s0[!is.na(s0$comp) & !is.na(s0$coeff),] 
s.sp <- sparseMatrix(i=as.numeric(s0$compound), j=as.numeric(s0$rxn), x=as.real(s0$coeff))
dimnames(s.sp) <- list(levels(s0$compound), levels(s0$rxn))
nmets <- nrow(s.sp); nrxns <- ncol(s.sp)
#code uses s as s.sp
s <- s.sp

##fba
fba.ub <- rep(10^3, nrxns); names(fba.ub) <- colnames(s.sp)
c.rxns <- c('SUCROSE-TRANS-RXN-L2R','CIT-TRANS-RXN-L2R', 'PROTON-TRANS-RXN-L2R')
fba.ub[vogel$rxn] <- 0
#fba.ub[c.rxns] <- 1000
fba.lp <- FBA(a=s.sp, sense='E')
v <- fba.lp$xopt; names(v) <- names(fba.obj)
#gapless
s.gapless <- s.sp
s.gapless[s.gapless!=0] <- s.gapless[s.gapless!=0]-10^-3
fba.lp <- FBA(a=s.gapless, sense='E')

##min flux adjust
#decision vars are [v r]
obj <- c(rep(0, nrxns), rep(1, nmets))
a <- cBind(s.sp, diag(nmets))
lb <- rep(0, nrxns+nmets)
ub <- rep(1000, nrxns+nmets)
names(obj) <- names(lb) <- names(ub) <- c(colnames(s.sp), rownames(s.sp))
lb['biomass'] <- 1
obj[as.character(s0$comp[s0$rxn=='biomass'])] <- 10^6
#run
cp <- Rcplex(cvec=obj, Amat=a, bvec=b, lb=lb, ub=ub, sense='G', objsense='min')
#examine r
x <- cp$xopt; names(x) <- names(obj)
r <- x[(nrxns+1):(nrxns+nmets)]
r[r>0]

##gapless
ss <- s.sp
ss[ss!=0] <- ss[ss!=0]-10^-6
fba.gl <- FBA(a=ss, fba.ub=fba.ub)
#goal prog
fba.goal <- FBA(a=ss, fba.obj.rxns=NULL, goal.rxns=unique(bm$rxn), fba.ub=fba.ub)
v <- fba.goal$xopt; names(v) <- colnames(ss); v[unique(bm$rxn)]
