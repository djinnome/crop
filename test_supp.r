##jmd
##10.7.11
##test_supp.r

#adds '[CCO-EXTRACELLULAR]' to add.supp
test.supp <- function(sp, ko.lst, ub=named.vec(10**3, names=colnames(sp)), supp.ub=10, annot, sense='E', quiet.fba=TRUE, add.extracellular=TRUE, ctrl=list(trace=0, method=1), ...){
    supps <- parse.supp(annot)
    #names
    nms <- intersect(intersect(names(ko.lst), rownames(annot)), names(supps))
    ko2 <- ko.lst[nms]; annot2 <- annot[nms,]; supps <- supps[nms]
    #stopifnot(all(names(ps) %in% rownames(annot.df)))
    mat <- matrix(NA, nrow=2, ncol=length(ko2), dimnames=list(c('nosupp', 'supp'), nms))
    supp.lst <- list()
    #loop thru ncu's
    for (i in 1:length(ko2)){
        if (annot2$ReplaceCPD[i]!=''){ rm.nuts <- annot2$ReplaceCPD[i] } else { rm.nuts <- NULL }
        if (annot2$ReplaceWithCPD[i]!=''){ add.nuts <- annot2$ReplaceWithCPD[i] } else { add.nuts <- NULL }
        ck0 <- check.new.media(sp, ub=ub, supp.ub=supp.ub, ko=ko2[[i]], rm.nuts=rm.nuts, add.nuts=add.nuts, ctrl=ctrl, quiet=TRUE, ...)
        if (ck0$status!=1){ cat(names(ko2)[i], 'has cplex stat', ck0$status, '\n') }
        mat[1, i] <- ck0$obj
                
        #loop thru potential nutrient conditions of ko i
        if (ck0$obj<10**-6){   
            ck.v <- rep(NA, length(supps[[i]]))
            for (j in 1:length(supps[[i]])){
                if (add.extracellular){ 
                    add.supp <- paste(supps[[i]][[j]], '[CCO-EXTRACELLULAR]', sep='') 
                } else { 
                    add.supp <- supps[[i]][[j]]
                }#end else
                add <- union(add.nuts, add.supp)
                
                if (all(add %in% rownames(sp))){
                    cnm <- check.new.media(sp, ko=ko2[[i]], ub=ub, supp.ub=supp.ub, rm.nuts=rm.nuts, add.nuts=add, sense=sense, quiet=quiet.fba, ctrl=ctrl, ...)
                    ck.v[j] <- cnm$obj
                    if (!(cnm$status %in% 1:2)){ cat(names(ko2)[i], add, 'has cplex stat', cnm$status, '\n') }
                } else { 
                    cat(names(ko2)[i], 'has supplements', add[!(add %in% rownames(sp))] ,'not in model \n') 
                }#end else
                
                supp.lst[[ nms[i] ]][[paste(add.supp, collapse=', ')]] <- ck.v[j]
                 
            }#end for j in supps
            mat[2, i] <- max(ck.v, na.rm=TRUE)
        }#end if ck0$obj
    }#end for i in ko's
    cat('\n Mutants that grow w/o supplements:', colnames(mat)[mat[1,]>10**-6], '\n')
    
    return(list(mat=mat, supp.lst=supp.lst))
}
