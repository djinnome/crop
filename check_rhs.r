##jmd
##6.14.2011
##check_rhs.r

library('Matrix')
library('Rcplex')
options(stringsAsFactors=FALSE)
setwd('/msc/neurospora/FBA/farm_data')
source('../farm/fba.r')

##made sparse format biomass_s.txt using make_biomass.r
##created sparse format transport rxns (import Vogel's media) in excel

##read
print("Checking growth on Vogels with Sv>=0")
bm <- read.delim('biomass_s.txt')
vogel <- read.delim('vogel_trans_s.txt')
s0 <- read.delim('s15.txt', as.is=FALSE, na='NIL')
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
fba.lp <- FBA(a=s.sp)
v <- fba.lp$xopt; names(v) <- names(fba.ub)

##get rhs
rhs <- (s.sp %*% v)[,1]
names(rhs) <- rownames(s.sp)
rhs.pos <- rhs[rhs>0]
