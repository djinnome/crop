	##jmd
##4.3.12
##xylose2etoh.r

#k is grams of neurospora / litre of soln
xylose2etoh <- function(s2, ub, suc.ub=0, k=3, xr.nad.ub=0.28/k, xylose.in.v=c(0.28, 0.84, 0.94, 1.5, 5.4)/k, o.thresh.v=c(0, 3.3, 5.1, 8.4, 12.6)/k, fba.obj='NGAM'){
    ub['SUCROSE-TRANS-RXN-L2R'] <- suc.ub
    lb <- rep(0, ncol(s2)); names(lb) <- colnames(s2); lb['NGAM'] <- 1.9
    #o.thresh.v <- c(seq(0, 0.1, by=0.001), seq(0.1, 1, by=0.1), seq(1, 10, by=0.1), seq(10, 100, by=10))
    #o.thresh.v=c(0, 10^(-4:3))
    
    #o.thresh.v <- 12.6+c(0, 3.3, 5.1, 8.4, 12.6)/k
    #xylitol.yield.v <- c(0, .16, .33, .06, 0)
    #xylitol.yield.v <- rep(0,5)
    
    mat <- data.frame(o.thresh.v, etoh=NA, etoh.yield=NA, o2=NA, xylose=NA, growth=NA, atp=NA)
    #rownames(mat) <- c(0, 3.3, 5.1, 8.4, 12.6)
    rownames(mat) <- o.thresh.v
    flux <- matrix(NA, ncol(s2), length(o.thresh.v), dimnames=list(colnames(s2), o.thresh.v))
    for (i in 1:length(o.thresh.v)){
        ub.tmp <- ub
        o.thresh <- o.thresh.v[i]
        ub.tmp['OXYGEN-MOLECULE-TRANS-RXN-L2R'] <- o.thresh
        ub.tmp['RXN2T-9-L2R[CCO-EXTRACELLULAR]'] <- xylose.in.v[i]
        ub.tmp['RXN2T-10-L2R'] <- xr.nad.ub
        #xylitol = 0.152 g/mmol; xylose = .150 g/mmol, so call it awash
        #lb['xylitol-out'] <- xylitol.yield.v[i]*ub.tmp['xylose-in']
        
        f <- FBA(s2, fba.ub=ub.tmp, fba.lb=lb, fba.obj.rxn=fba.obj, quiet=TRUE, min2=FALSE, control=list(trace=0))
        flux[,i] <- f$xopt
        if (f$stat==1){
            #etoh is 46 g/mol; xylose is 150 g/mol
            mat[i, 'etoh'] <- f$x['ETOH-TRANS-RXN-R2L']
            mat[i, 'etoh.yield'] <- 46*f$x['ETOH-TRANS-RXN-R2L'] / (150*f$xopt['xylose-in'])
            mat[i, 'o2'] <- f$xopt['TRANS-RXN2T-231-L2R[CCO-EXTRACELLULAR]']
            mat[i, 'xylose'] <- f$xopt['RXN2T-9-L2R[CCO-EXTRACELLULAR]']
            #mat[i, 'xylitol'] <- f$xopt['xylitol-out']
            mat[i, 'growth'] <- f$xopt['biomass']
            mat[i, 'atp'] <- f$x['NGAM']
        }
        #cat(o.thresh, 'gives etoh', f$x['etoh-out'], '\n')
        #cat('flux thru etoh is: \n'); print(fom <- flux.of.met(s2, f$x, 'ETOH'))
    }
    mat <- mat[order(mat$etoh.yield),]
    #print(conversion <- mat$etoh/(1.67*mat$xylose))
    return(list(output=signif(mat, digits=3), flux=flux))
}
