##jmd
##11.13.12
##explain_SLs_fcns.r

## loop thru all gene pairs from pwys.sl involving gene.name, and run FBA w/ export of met.export
#takes gpr from larger environment
testGeneOfSLwExport <- function(sp, gene.name, met.export, pwys.sl){    
    # get subset pwys.tmp
    pwys.tmp <- pwys.sl[which(pwys.sl[,'name1'] %in% gene.name | pwys.sl[,'name2'] %in% gene.name),]
    print(pwys.tmp[,1:4])
    # set up constraint senses
    se=named.vec('E', rownames(sp))
    se[met.export] <- 'G'
    # loop through gene pairs
    g <- named.vec(NA, rownames(pwys.tmp))
    for (i in 1:nrow(pwys.tmp)){
        ko.rxns <- ncu2rxn(c(pwys.tmp$Gene1[i], pwys.tmp$Gene2[i]), gpr=gpr)
        f <- FBA(sp, ko=ko.rxns, fba.ub=fva.ub, se=se, control=list(trace=FALSE))
        g[i] <- f$obj
    }
    return(g)
}

#take gpr from env
wrapMinExport <- function(sp, ncu.v, fva.ub=named.vec(1000, colnames(sp))){
    pen.mets <- character(0)
    w <- named.vec(1, rownames(s.al))
    ko.rxns <- ncu2rxn(ncu.v, gpr=gpr)
    f <- FBA(s.al, ko=ko.rxns, fba.ub=fva.ub, se='E')
    if (f$obj>0.01){
        print('no need for min.export')
    } else {
    cond=FALSE
    me=named.vec(0, 'first')
    s2 <- s.al[,setdiff(colnames(s.al), ko.rxns)]

        while (length(me)==1 && !cond){
            me <- min.export(s2, ub=fva.ub[colnames(s2)], v.min=0.1, w=w)
            w[names(me)] <- 10**3

        if (length(pen.mets)>0){
        cond=names(me) %in% pen.mets
        }
        if (length(me)==1) pen.mets <- c(pen.mets, names(me))
        }#end while
    }
    return(pen.mets)
}

getAssocGenes <- function(pwys.sl, gene.name){
    gene2.v <- unlist(apply(pwys.sl, 1, FUN=function(x){ 
        if (x['name1'] %in% gene.name | x['name2'] %in% gene.name)
            setdiff(x[c('name1', 'name2')], gene.name) 
        }
    ))
    return(gene2.v)
}
