##jmd
##6.5.12
##multiGeneKO.r
##KO 2 genes & run FBA

multiGeneKO <- function(a, ncu.mat, gpr, max.fba=1, ...){
    ncu.mat <- as.matrix(ncu.mat)
    cat('There are', ncol(ncu.mat), 'columns. 2 imples double KO.\n')
    obj.v <- apply(ncu.mat, 1, FUN=function(x){
        ko.rxns.tmp <- ncu2rxn(ncu=x, gpr=gpr)
        if (length(ko.rxns.tmp)>=1){
            fba.tmp <- fba.na(a=a, ko=ko.rxns.tmp, ...)
            ret <- signif(fba.tmp$obj, 3)
        } else {
            ret  <- max.fba
        }#end else
    }) #end apply
    names(obj.v) <- rownames(ncu.mat)
    return(obj.v)
}
