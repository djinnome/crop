##jmd
##2.13.11
##make_rxn_annot.r
##takes in S and annot files & returns rxn.annot
#use s15 since don't need balancing for fba w/ v_bm<=1

make.rxn.annot <- function(sp, univ.ec, org.ec, biomass.rxns, pwy.rxns, suffix.regexp='\\[CCO-.+\\]', known.mets.conc=NULL, gpr, gpr.rxn.col='RXN', gpr.gene.col=3, gpr.spont.gene='s0001', org.name='nc',
eqn.col='eqn', kegg.col='Kegg', rxn.col='rxn', ecpred.file='/msc/neurospora/FBA/Neurospora/NC10_CALLGENES_FINAL_2.ecpred', rxns.imb.file='rxns_imbalance.csv', write.file='rxn_annot.txt'){

    bal.mat <- read.csv(rxns.imb.file, row.names=1)
    #get eficaz probs
    prob.mat <- ef.probs(ecpred.file, low.conf.prob=0.3)
    #write.csv(round(ef.prob.v, 2), 'ec_probs.csv')
    
    ##clean meta.ec kegg/eqn columns
    univ.ec[,eqn.col] <- gsub(' ', '', univ.ec[,eqn.col])
    univ.ec[univ.ec[,eqn.col]=='', eqn.col] <- NA
    univ.ec[,kegg.col] <- gsub(' ', '', univ.ec[,kegg.col])
    univ.ec[univ.ec[,kegg.col]=='', kegg.col] <- NA
    univ.ec[is.na(univ.ec[,kegg.col]), kegg.col] <- univ.ec[is.na(univ.ec[,kegg.col]), eqn.col]
    univ.ec$Reason <- univ.ec$Pathways <- NA
    
    #combine
    #1.30.12: prioritize nc.ec
    rxn.ec <- merge.ec(univ=univ.ec, org=org.ec, rxn.col=rxn.col, priority='org')
    #match
    rxn.ec.match <- rxn.ec[colnames(sp),]
    rxn.ec.match$rxn <- rownames(rxn.ec.match) <- colnames(sp)
    
    ##annotate nc status
    rxn.ec.match[[org.name]] <- as.numeric(rxn.ec.match$rxn %in% c(biomass.rxns, org.ec[,rxn.col]))
    #nc.pwy status, independently of cellular compartment
    rxns.nosuffix <- gsub(suffix.regexp, '', rxn.ec.match[,rxn.col])
    #count the number of pwy's a rxn is in using pwy.rxns column
    rxn.ec.match$pwy <- apply(as.matrix(rxns.nosuffix), 1, FUN=function(x){ sum(x==pwy.rxns, na.rm=TRUE) })
    
    ##parse ec numbers
    rxn.ec.match$EC <- gsub('EC# ', '', rxn.ec.match$EC)
    #extract ec from rxn where available, else ec=NA
    #some rxns have ec in rxn name but not in ec column
    hidden.ec <- setdiff(grep('.+\\..+\\..+\\..+', rxn.ec.match$rxn), grep('EC#', rxn.ec.match$EC))
    rxn.ec.match$EC[hidden.ec] <- gsub('-RXN.+', '', rxn.ec.match$rxn[hidden.ec])
    #some have rxn name in ec column
    rxn.ec.match$EC[-grep('.+\\..+\\..+', rxn.ec.match$EC)] <- NA
    
    ##match probs to annotations
    rxns.annot <- merge(x=rxn.ec.match, y=prob.mat, all.x=TRUE)
    rownames(rxns.annot) <- rxns.annot[,rxn.col]
    rxns.annot <- rxns.annot[colnames(sp),]
    stopifnot(all(rownames(rxns.annot)==rownames(bal.mat)), all(rownames(rxns.annot)==colnames(sp)))
    
    #unbalanced rxns
    #bal.mat <- bal.mat
    #catch rxns that generate mass
    rxns.annot$h.unb <- bal.mat[,'H']
    rxns.annot$unb.add <- as.numeric(rowSums(bal.mat[,c('C', 'N', 'O', 'P', 'S')]>0))
    rxns.annot$unb.rm <- as.numeric(rowSums(bal.mat[,c('C', 'N', 'O', 'P', 'S')]<0))
    
    #rev rxns
    bare.rxns <- gsub('-L2R$|-R2L$', '', rxns.annot[,rxn.col])
    rxns.annot$rev <- as.numeric(bare.rxns %in% bare.rxns[duplicated(bare.rxns)])
    
    ##gpr
    rxns.annot$gpr <- as.numeric(rxns.annot[,rxn.col] %in% gpr[,gpr.rxn.col]) #from 3-column GPR
    rxns.annot$spont <- as.numeric(rxns.annot[,rxn.col] %in% gpr[gpr[,gpr.gene.col] %in% gpr.spont.gene, gpr.rxn.col])
    
    ##thermo probs
    rxns.annot$dg0 <- get.dg0(s=sp)
    rxns.annot$thermo.prob <- dg2prob(s=sp, dg0.rxns=rxns.annot$dg0, known.mets.conc=known.mets.conc)
    #rxns$thermo.prob <- NA
    
    ##write
    #put eqn as last column
    rxns.out <- data.frame(rxns.annot[,colnames(rxns.annot)!='eqn'], eqn=rxns.annot[,colnames(rxns.annot)=='eqn'])
    rownames(rxns.out) <- rxns.out[,rxn.col]    
    if (!is.null(write.file)){ write.table(rxns.out, write.file, sep='\t', quote=FALSE) }
    return(rxns.out)
}

###check#########################################################################################
