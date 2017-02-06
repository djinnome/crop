##jmd
##2.2.12
##call_farm.r

##to do: need to add known mets, consider reverse of rxn direction of irrev rxns, based on thermo

source('/msc/neurospora/FBA/farm/farm_init.r')

fva.ub2 <- fva.ub; fva.ub2["AMMONIUM-TRANS-RXN-L2R"] <- 0
f.ub.mat <- as.matrix(cbind(fva.ub, fva.ub2))
nconds <- ncol(f.ub.mat)

#########################################################################################################################
####get beta0
#########################################################################################################################
##1 condn
x0 <- get.x0.md.multi(sp=sp.f, obj.coeffs=p-0.9, biomass='biomass', bm.ub=0.3, lb.req=0.1, check.input=TRUE, se='E', f.ub.mat=f.ub.mat, mipstart=NULL,
md=FALSE, known.rxns=unlist(known.rxns), nondilute.regexp=nondilute.regexp, ctrl=list(tilim=300, numericalemphasis=TRUE), connect.constr.mat=cut.met.mat)

bv <- x0$mat[,'beta.v']; names(bv) <- gsub('cond1_', '', names(bv))
beta0 <- x2beta(sp=sp.f, x0=bv, obj=p-0.9, f.ub=fva.ub)

##add md constraints
#don't need met.ub since this analysis has its own large bounds and should apply to *any* media

##all
mipstart <- c(numeric(ncol(sp.f)*nconds), as.numeric(colnames(sp.f) %in% rownames(beta0)))
x0 <- get.x0(sp=sp.f, obj.coeffs=p-0.9, biomass='biomass', bm.ub=0.3, lb.req=0.1, check.input=TRUE, se='E', f.ub.mat=f.ub.mat, mipstart=mipstart, known.rxns=known.rxns,
md=TRUE, nondilute.regexp=nondilute.regexp, ctrl=list(tilim=3600, numericalemphasis=TRUE), connect.constr.mat=cut.met.mat, ko.rxns.lst=list('OXYGEN-MOLECULE-TRANS-RXN-L2R'))

bv <- x0$mat[,'beta.v']; names(bv) <- gsub('cond._', '', names(bv))
beta0 <- x2beta(sp=sp.f, x0=bv, obj=p-0.9, f.ub=fva.ub)

#########################################################################################################################
####check
#########################################################################################################################
beta1 <- rownames(beta0)
mean(beta1 %in% model.rxns)
unlist(nut.rxns) %in% beta1
ff <- FBA(sp.f[,beta1], fba.ub=fva.ub[beta1])

check.ko(sp.f[,beta1], ko.lst=nut.rxns, ub=fva.ub[beta1])
