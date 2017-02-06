##jmd
##2.18.11
##make_dual_mat.r
##for no-growth conditions

library('Rcplex')

setwd('/home/radon00/jdreyf/farm')
#setwd('Z:/meta_recon')

##read
#all irrev
s <- read.csv('test_s_irrev.csv', row.names=1, as.is=FALSE)
s <- as.matrix(s)+1-1 #make 'real' matrix
nmets <- nrow(s); nrxns <- ncol(s)
rxn.rev <- rep(FALSE, nrxns)
rxn.rev[9:10] <- TRUE
rev.ind <- which(rxn.rev)

##ko
ko.ind <- 4

##primal vars
primal.obj <- primal.lb <- rep(0, ncol(s))
primal.ub <- rep(1000, ncol(s))
names(primal.obj) <- names(primal.lb) <- names(primal.ub) <- colnames(s)
#set values
lim.nut.ind <- c(1,5,10)
lim.nut.max <- c(1,0,1)
primal.ub[lim.nut.ind] <- lim.nut.max
primal.ub[ko.ind] <- 0
primal.obj['biomass'] <- 1
#lp
(p <- Rcplex(cvec=primal.obj, Amat=s, bvec=rep(0, nrow(s)), lb=primal.lb, ub=primal.ub, objsense='max', sense='E'))

##dual
#decision vars are [mu, l=lambda] of length nmets+nrxns
sign.lambda <- rep(1, nrxns)
#nuts & ko rxns face *upper* bounds, opposite of other fluxes
sign.lambda[lim.nut.ind] <- -1
a <- cbind(t(s), diag(sign.lambda))
rhs <- -primal.obj
lb <- rep(c(-Inf, 0), times=c(nmets, nrxns))
ub <- rep(c(Inf, Inf), times=c(nmets, nrxns))
names(lb) <- names(ub) <- c(paste('mu',1:nmets,sep=''), colnames(s))
#duals of ko'd are free
lb[nmets+ko.ind] <- -Inf
#treating rev rxns as 2 independent irrev rxns works, but treating them as reversible rxns 
#by ignoring the two v>=0 constraints does not
#b/c one rxn may want to go below 0 even if other is at its max (imposed by nutrients)
#ub[nmets+setdiff(rev.ind, union(ko.ind, lim.nut.ind))] <- 0
dual.obj <- rep(0, nmets+nrxns)
dual.obj[nmets+lim.nut.ind] <- lim.nut.max
#lp
(d <- Rcplex(cvec=dual.obj, Amat=a, bvec=rhs, sense='E', lb=lb, ub=ub))
(mu <- d$x[1:nmets])
lam <- d$x[(nmets+1):(nmets+nrxns)]
names(lam) <- colnames(s)
lam
all.equal(d$obj,p$obj)
            
###check#################################################################################################################
##checked both primal & dual using different conditions & diff ko, & both had same obj value
