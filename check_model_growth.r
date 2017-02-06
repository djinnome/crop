##jmd
##8.22.2011
##check_model_growth.r

source('/msc/neurospora/FBA/farm/farm_config.r')

##make rxns
rxns <- data.frame(matrix(0, nrow=ncol(s.sp), ncol=2, dimnames=list(colnames(s.sp), c('nc', 'prob'))))
rxns[rownames(rxns) %in% c(nc$rxn, bm$rxn, bm.goals$rxn, vogel$rxn), 'nc'] <- 1
#check
print('Check that annotations match matrix')
all(colnames(s.sp)==rownames(rxns))

##subset
model.rxns <- read.csv('model_rxns.csv')[,1]
ss <- intersect(model.rxns, colnames(sg))
s.al <- sg[,ss]
s.al <- s.al[rowSums(abs(s.al))>0,]
rxns <- rxns[ss,]
nmets <- nrow(s.al); nrxns <- ncol(s.al)

##fba
fva.ub <- rep(10**5, nrxns); names(fva.ub) <- colnames(s.al)
fva.ub[intersect(vogel$rxn, colnames(s.al))] <- 1000
#get glucose uptake of ~18 in yeast, which corresponds to sucrose uptake of ~9, from snitkin et al (http://www.jsbi.org/pdfs/journal1/IBSB08/IBSB08011.pdf)
#fva.ub['SUCROSE-TRANS-RXN-L2R'] <- 9

print("Checking growth of model on Vogels with Sv>=0 and gapless constraints")
fba.lp <- FBA(a=s.al, control=list(trace=0), fba.ub=fva.ub, sense='G', quiet=FALSE)

if (fba.lp$obj>10**-6){
    print("Checking growth of model on Vogels with Sv=0 and gapless constraints")
    fba2 <- FBA(a=s.al, control=list(trace=0),  fba.ub=fva.ub, sense='E', quiet=FALSE)
}

##min flux adjust
if (fba.lp$obj<10^-6){
    print("Running min metabolite adjustment to analyze failure")
    min.met.adj(s.al)
}
