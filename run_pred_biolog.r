##jmd
##3.4.12
##run_pred_biolog.r

#sp=s.al; default.nuts=vogel
#add.trans=TRUE; exclude.trans.mets='CIT'; plate.names='PM1';
#extracellular=TRUE; write.out.tsv=NULL; na.is.nogrowth=TRUE;
#biol.nuts.file='biolog.nutrients'; biol.wells.file='biolog_wells.tsv'; ko.rxns=NULL
run.pred.biolog <- function(sp, default.nuts, add.trans.mets=NULL, plate.names=paste('PM', 1:4, sep=''),  exclude.trans.mets=NULL, na.is.nogrowth=TRUE, 
extracellular=TRUE, write.out.tsv='biolog_predict', biol.nuts.file='biolog.nutrients', biol.wells.file='biolog_wells.tsv', fba.sense='E', eps=10**-6,
ko.rxns=NULL, supp.ub=10){  
    bw <- read.delim(biol.wells.file)
    bw <- bw[rowSums(bw!='')>0,]
    #keep negative controls but get rid of other ReplaceWithCPD=''
    #bw <- bw[union(which(bw$ReplaceWithCPD!=''), grep('negative control', bw$MediumName)),]
    
    bw$PlateID <- gsub('Biolog | - .+Sources$', '', bw$PlateID)
    bw <- bw[order(bw$PlateID, bw$WellId),]
    bw <- bw[bw$PlateID %in% plate.names,]
    
    #strsplit(x='', split=',') gives 'character(0)'
    model.mets <- unique(gsub('\\[CCO-.+\\]$', '', rownames(sp)))
    bw <- bw[sapply(strsplit(bw$ReplaceWithCPD, ','), FUN=function(x){ any(x %in% model.mets | length(x)==0) }),]
    
    #add [CCO-EXTRACELLULAR]?
    if (extracellular){ 
        bw$ReplaceWithCPD <- apply(bw, MARGIN=1, FUN=function(x){ 
            if (x['ReplaceWithCPD']!=''){ 
                paste(paste(unlist(strsplit(x=x['ReplaceWithCPD'], split=',')), '[CCO-EXTRACELLULAR]', sep=''), collapse=',')
            } else { x['ReplaceWithCPD'] }#end else
        })
    }#end if extracellular
    
    #?
    if (length(add.trans.mets)>0){
         bw$ReplaceWithCPD <- apply(bw, MARGIN=1, FUN=function(x){ 
            y <- x['ReplaceWithCPD']
            if (y!=''){
                add.nuts.v <- gsub('\\[CCO-.+\\]$', '', unlist(strsplit(x=y, split=',')))
                if (any(add.nuts.v %in% add.trans.mets)){
                    y[add.nuts.v %in% add.trans.mets] <- add.nuts.v[add.nuts.v %in% add.trans.mets]
                    y <- paste(y, collapse=',')
                } else { y } #end if any
            } else { y } #end if is empty
        }) #end apply fcn
    } #end if add.trans.mets
        
    bn <- read.delim(biol.nuts.file)
    pb <- predict.biolog(sp=sp, default.nuts=default.nuts, extracellular=extracellular, write.out.tsv=write.out.tsv, bn=bn, bw=bw, fba.sense=fba.sense, ko.rxns=NULL, supp.ub=supp.ub, eps=eps)
    if (na.is.nogrowth){ pb$pred[is.na(pb$pred)] <- 0 }
    return(pb)
}
