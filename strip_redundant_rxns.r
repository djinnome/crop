##jmd
##9.17.11
##strip_redundant_rxns.r

strip.redundant.rxns <- function(sp, rxns.annot){
    stopifnot(rownames(rxns.annot)==colnames(sp))
    #condense columns into strings, normalizing smallest |coeff|=1
    rxns.v <- apply(sp, 2, FUN=function(x){ x <- x[x!=0]; paste(names(x), x/min(abs(x)), sep=':', collapse=';') })
    #order rxns by: instantiated, prob, ...
    rxns.v.o <- rxns.v[order(!(names(rxns.v) %in% grep('//', names(rxns.v), value=TRUE)), rxns.annot$nc, rxns.annot$gpr, rxns.annot$nc.pwy, rxns.annot$prob, rxns.annot$rev, -nchar(names(rxns.v)), decreasing=TRUE)]
    #delete dups
    rxns.out <- rxns.v.o[!duplicated(rxns.v.o)]
}

get.redundant.cytosolic.rxns <- function(rxns, model){
    #intersect ( model.rxns not in nccyc ) with ( location-free description of nccyc rxns )
    (ret <- intersect(intersect(model, rownames(rxns)[rxns$nc==0]), gsub('\\[CCO-.+\\]', '', rownames(rxns)[rxns$nc==1])))
    #check by looking at rxns in ret & their dups in nccyc
    print(rownames(rxns)[gsub('\\[CCO-.+\\]', '', rownames(rxns)) %in% ret])
    cat('Done check \n')
    return(ret)
}

#create strings of rxns, replace nadp by nad & nadph by nadh in instantiated rxns whose generic rxn has nadp-or-nop | Acceptor | Donor-H2
strip.nadp <- function(sp, rxns.annot, model.rxns){
    nads <- c('NAD', 'NADP', 'NADH', 'NADPH')
    
    #want to put nads at end, so can replace among them in strings w/o worry about order of mets
    sp <- sp[c(setdiff(rownames(sp), nads), nads), c(intersect(model.rxns, colnames(sp)), setdiff(colnames(sp), model.rxns))]
    #condense columns into strings, normalizing smallest |coeff|=1
    rxns.v <- apply(sp, 2, FUN=function(x){ x <- x[x!=0]; paste(names(x), x, sep=':', collapse=';') })
    
    #get generics, their frames & rxns w/ nad(p), which shouldn't be in model.rxns
    generics <- grep('(^|;)(NADP-OR-NOP|Acceptor|Donor-H2):', rxns.v, value=TRUE)
    gen.frames <- rxns.annot[names(generics), 'FRAME.ID']
    instant.rxns.of.generics <- rxns.v[intersect(grep('//', model.rxns, value=TRUE), rownames(rxns)[rxns$FRAME.ID %in% gen.frames])]
    
    #get instantiated model rxns w/ nadp(h), which is always at end of rxns.v string
    inst.nadp <- rxns.v[ names(instant.rxns.of.generics)[grep(';NADP(H|):', instant.rxns.of.generics)] ]
    #generate strings using nad(h) instead of nadp(h), to see if they're in rxns.v
    inst.nad <- gsub(pattern='NADP', replacement='NAD', gsub(pattern='NADPH', replacement='NADH', inst.nadp))
    #new set of rxns to replace nadp(h) ones
    rxns.annot.v <- rxns.v[ rownames(rxns.annot)[rxns.annot$nc==1] ]
    rxns.nad <- names(rxns.annot.v)[match(inst.nad, rxns.annot.v)]; names(rxns.nad) <- names(inst.nad)

    return(rxns.nad)
}

#use:
#model.rxns <- read.csv('model_rxns.csv')[,1]
#sn <- strip.nadp(sp=s.sp, rxns.annot=rxns, model.rxns=colnames(s.al))
#model.rxns <- union(setdiff(model.rxns, names(sn)), sn)

##replace nadp w/ nad using classes
#strip.nadp <- function(sp, nc.pwy, pwy.classes){
#    nads <- c('NAD','NADP','NADH','NADPH')    
#    sn <- sp[!(rownames(sp) %in% nads), colSums(sp[nads,]!=0)==2]
#    sn <- as.matrix(sn[,colSums(abs(sn))>0])
#    dups <- sn[,duplicated(t(sn))]
#    colnames(dups)
#}
