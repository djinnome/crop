##jmd
##2.9.12
##x2beta.r

#vectors must be named & match colnames(sp)
x2beta <- function(sp, x0, obj.coeff, thresh=c(0.9, 0.5, 10**-(6:9), 0), biomass='biomass', f.ub=array(10**3, dim=ncol(sp), dimnames=list(colnames(sp)))){
    stopifnot(length(x0)==ncol(sp), names(x0)==colnames(sp), names(obj.coeff)==colnames(sp), names(f.ub)==colnames(sp))
    thresh <- thresh[order(-thresh)]
    
    #check growth of model.rxns & new rxns combined
    for (thr in thresh){
        rxns1 <- colnames(sp)[x0>=thr]
        if (biomass %in% rxns1){
            cat('Testing if new rxn set gives FBA growth for Sv>=0 w/ threshold', thr, ': ')
            fba.svg0 <- FBA(sp[,rxns1], fba.ub=f.ub[rxns1], control=list(trace=0), sense='G')
            
            cat('Testing if new rxn set gives FBA growth for Sv=0 w/ threshold', thr, ': ')
            fba <- FBA(sp[,rxns1], fba.ub=f.ub[rxns1], control=list(trace=0), sense='E')
            cat('\n')
            
            if (fba$obj>10**-6){
                ret <- matrix(cbind(x0[rxns1], obj.coeff[rxns1], fba$xopt), ncol=3, dimnames=list(rxns1, c('beta', 'obj.coeff', 'v')))
                (ret.o <- ret[order(-ret[,1]),])
                break
            }
        }#end if biomass
    }#end for
    
    return(ret.o)
}
