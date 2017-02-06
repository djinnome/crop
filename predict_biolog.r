##jmd
##1.5.12
##predict_biolog.r

#NOTE: changeNuts() automatically treats rm.nuts as extracellular BUT makes no assumptions about compartment of add.nuts.
#this pulls object, e.g. model.rxns, from main environment
predict.biolog <- function(sp, default.nuts, bn, bw, extracellular=FALSE, write.out.tsv='biolog_predict', fba.sense='E', ko.rxns=NULL, supp.ub=10, eps=10**-6){  
    ##change default medium
    rm.nuts <- gsub('\\[CCO-EXTRACELLULAR\\]$', '', setdiff(default.nuts$comp, bn$comp))
    #cat('Checking growth on Biolog default medium:', '\n')
    s2 <- changeNuts(sp, rm.nuts=rm.nuts, add.nuts=setdiff(bn$comp, default.nuts$comp))
    nmets <- nrow(s2); nrxns <- ncol(s2)
    #ub
    fva.ub0 <- named.vec(10**3, colnames(s2))
    fva.ub0[intersect(c('PYRUVATE-TRANS-RXN-L2R', 'AMMONIUM-TRANS-RXN-L2R'), names(fva.ub0))] <- c(25, 50)
    #f <- FBA(s2, sense='E', fba.ub=fva.ub0, ko.rxns=ko.rxns, control=list(trace=0))
    
    ##annotate
    plate <- gsub('Biolog | - .+Sources$', '', bw$Plate)
    pred <- quant.pred <- rep(NA, nrow(bw)); names(pred) <- paste(plate, bw$ReplaceWithCPD, sep='.')
    
    ##loop thru wells & predict
    for (i in 1:nrow(bw)){
        rm.nuts <- unlist(strsplit(bw$ReplaceCPD[i], ','))
        add.nuts <- intersect(unlist(strsplit(bw$ReplaceWithCPD[i], ',')), rownames(s2))
        if (length(add.nuts)>0){
            fba.tmp <- check.new.media(sp=s2, ub=fva.ub0, supp.ub=supp.ub, rm.nuts=rm.nuts, add.nuts=add.nuts, sense=fba.sense)
            stopifnot(fba.tmp$stat %in% 1:2)
            quant.pred[i] <- fba.tmp$obj
        }
    }#end for
    quant.pred[quant.pred<0] <- 0
    #leave preds between 0 and 10**-6 as is
    pred[!is.na(quant.pred) & quant.pred>eps] <- 1
    pred[!is.na(quant.pred) & quant.pred<0] <- 0
    
    ##return
    df.out <- data.frame(plate=plate, well=bw$Well, medium=bw$MediumName, cpd=bw$ReplaceWithCPD, pred=pred, quant.pred=quant.pred)
    rownames(df.out) <- paste(df.out$plate, df.out$well, sep='.')
    df.out <- df.out[order(df.out$plate, df.out$well),]
    
    if (!is.null(write.out.tsv) && !is.na(write.out.tsv) && write.out.tsv!=''){
        if (extracellular){
            write.table(df.out, paste(write.out.tsv, '.tsv', sep=''), sep='\t', row.names=FALSE)
        } else {
            write.table(df.out, paste(write.out.tsv, '_intracellular.tsv', sep=''), sep='\t', row.names=FALSE)
        }#end else
    }#end if write
    
    return(df.out)
}#end fcn
