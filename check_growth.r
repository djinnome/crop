##jmd
##5.23.2011
##check_growth.r

source('/msc/neurospora/FBA/farm/farm_config.r')

##fba
fba.ub <- rep(10^3, nrxns); names(fba.ub) <- colnames(s.sp)
print("Checking growth on Vogels with Sv=0")
fba.lp <- FBA(a=s.sp, control=list(trace=0))
v <- fba.lp$xopt; names(v) <- names(fba.ub)

#gapless
print("Checking growth on Vogels with gapless constraints")
fba.lp <- FBA(a=sg, control=list(trace=0))
v <- fba.lp$xopt; names(v) <- colnames(s.sp)

##min flux adjust
if (fba.lp$obj<10^-6){
    print("Running min flux adjustment to analyze failure")
    min.flux.adj(sg)
}
