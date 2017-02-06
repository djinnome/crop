##jmd
##2.8.11
##check_ko.r

#... to FBA
check.ko <- function(s.test, rxn.names=colnames(s.test), ko.lst, sense='E', ub=rep(10**3, ncol(s.test)), obs=0, print.v=FALSE, annot.df=NULL, quiet.fba=TRUE, ctrl=list(trace=0, method=1), ...){
    if (!is.null(annot.df)){ 
        stopifnot(names(ko.lst) %in% rownames(annot.df)) 
        if (sum(annot.df$ReplaceWithCPD!='')>0){
            if(!all(paste(annot.df$ReplaceWithCPD[annot.df$ReplaceWithCPD!=''], '[CCO-EXTRACELLULAR]', sep='') %in% rownames(s.test))){ 
                stop('Not all replace CPDs in S') 
            }
        }
        annot.df <- annot.df[names(ko.lst),]
    }
    ko <- data.frame(name=names(ko.lst), obs=obs, pred=NA, all.inc=FALSE)
    
    ##set-up LP
    nmets <- nrow(s.test)
    nrxns <- ncol(s.test)
    names(ub) <- rxn.names
    #sense
    se <- character(nrow(s.test)); names(se) <- rownames(s.test)
    se[1:nrow(s.test)] <- sense
    
    ##loop thru ko
    for (i in 1:nrow(ko)){
        se.tmp <- se
        rxns.ko.tmp <- ko.lst[[i]]
        #address 'unable to use' phenotypes by 1. taking away ReplaceCPD and 2. adding ReplaceWithCPD.
        if (!is.null(annot.df)){ 
            if (!is.na(annot.df$ReplaceCPD[i]) && annot.df$ReplaceCPD[i]!=''){ 
                rxns.ko.tmp <- union(rxns.ko.tmp, paste(annot.df$ReplaceCPD[i], 'TRANS-RXN-L2R', sep='-')) 
            }
            if (annot.df$ReplaceWithCPD[i]!=''){ 
                se.tmp[paste(annot.df$ReplaceWithCPD[i], '[CCO-EXTRACELLULAR]', sep='')] <- 'L'
            }
        }#end if !is.null(annot.df)
        int <- intersect(rxns.ko.tmp, rxn.names)
        #if KOs
        if (length(int)>=1){
            if (all(rxns.ko.tmp %in% rxn.names)) ko[i, 'all.inc'] <- TRUE
            ub.ko <- ub
            #KO these rxns by setting ub to 0 (all rxns >= 0 since 2 rxns per reversible one)
            ub.ko[int] <- 0
            #fba
            #cp.tmp <- Rcplex(cvec=obj, Amat=s.test, bvec=rep(0, nrow(s.test)), lb=rep(0, nrxns), ub=ub.ko, objsense='max', sense=sense, control=list(trace=0))
            cp.tmp <- FBA(a=s.test, fba.ub=ub.ko, sense=se.tmp, control=ctrl, quiet=quiet.fba, ...)
            #check if growth
            if (cp.tmp$status==1){ ko[i, 'pred'] <- cp.tmp$obj } else { print(paste(ko$name[i], 'on row', i, 'has cplex status', cp.tmp$status)) }
            #print fluxes if not consistent
            pred.binary <- ko[i, 'pred']>=10^-6
            if (pred.binary!=obs[i] & print.v){ print.flux(v=cp.tmp$xopt, out.file=paste(ko$name[i], 'fluxes.tsv', sep='_')) }#end if
        }#end if length
    }#end for
    return(ko)
}

#v=cp.tmp$xopt
#names(v)=colnames(s.test)
#met <- 'XANTHINE'
#cbind(name=jdc(s.test, r=met), flux=v[names(jdc(s.test, r=met))])
