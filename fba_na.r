##jmd
##6.5.12
##fba_na.r
##run fba, but try diff methods if is.na(fba$obj)

fba.na <- function(fba.lb=NULL, ngam.val=0, control=list(trace=0), ...){
    fba.tmp <- FBA(fba.lb=fba.lb, ngam.val=ngam.val, control=control, ...)
    if (is.na(fba.tmp$obj) & (is.null(fba.lb)|all(fba.lb==0)) & ngam.val==0){
        rc.meth <- 0
        while (is.na(fba.tmp$obj) & rc.meth<=4){
            #control <- control[-which(names(control)=='method')]
            fba.tmp <- FBA(fba.lb=fba.lb, ngam.val=ngam.val, control=list(method=rc.meth, numericalemphasis=1, trace=0), ...)
            rc.meth <- rc.meth+1
        }#end while
    }#end if
    return(fba.tmp)
}#end fcn
