##jmd
##2.12.12
##FBAmodel-class.r

# http://stackoverflow.com/questions/10961842/how-to-define-the-subset-operators-for-a-s4-class

library('Biobase')
library('Matrix')

##names.ann = model annotation
#rxns needs to have cols: ub, lb, obj.c
#rxns & mets need to match S
#if gpr is non-empty, it needs to have cols: rxn, enzyme, gene
setClass(Class='FBAmodel', representation=list(S='Matrix', rxnData='AnnotatedDataFrame', metData='AnnotatedDataFrame', geneData='AnnotatedDataFrame', enzymeData='AnnotatedDataFrame', gpr='AnnotatedDataFrame'))

#setMethod(f='[', 

##accessor fcns
#don't want fcns w. common names like: rxns, mets, genes, since they'll class w/ names i already use; but don't want to have to call pData, either
#these can be used to extract but not to replace?
geneData <- function(fbam){
    stopifnot(class(fbam)=='FBAmodel')
    return(pData(fbam@geneData))
}

metData <- function(fbam){
    stopifnot(class(fbam)=='FBAmodel')
    return(pData(fbam@metData))
}

rxnData <- function(fbam){
    stopifnot(class(fbam)=='FBAmodel')
    return(pData(fbam@rxnData))
}

#names
rxnNames <- function(fbam){
    stopifnot(class(fbam)=='FBAmodel')
    return(colnames(fbam@S))
}

metNames <- function(fbam){
    stopifnot(class(fbam)=='FBAmodel')
    return(rownames(fbam@S))
}

geneNames <- function(fbam){
    stopifnot(class(fbam)=='FBAmodel')
    return(rownames(pData(fbam@geneData)))
} 

enzymeNames <- function(fbam){
    stopifnot(class(fbam)=='FBAmodel')
    return(rownames(pData(fbam@enzymeData)))
} 

#rxnData cols
ub.fbam <- function(fbam){ return(rxnData(fbam)$ub) } 
lb.fbam <- function(fbam){ return(rxnData(fbam)$lb) }
obj.fbam <- function(fbam){ return(rxnData(fbam)$obj) }

##fcns
#this doesn't rm a rxn when one of its mets are removed!
ss.fbam <- function(fbam, c=NULL, r=NULL){
    stopifnot(class(fbam)=='FBAmodel')
    sp <- fbam@S
    #can't use ss.mat, since this fcn also uses 'matrix' class
    #fbam@S <- ss.mat(sp, c, r)
    pData(fbam@rxnData) <- rxnData(fbam)[colnames(fbam@S),]
    pData(fbam@metData) <- metData(fbam)[rownames(fbam@S),]
    return(fbam)
}

##example
#(m <- Matrix(c(0,2,0,1,0,-1), nrow=2, dimnames=list(c('r1', 'r2'), letters[1:3]), byrow=TRUE))
#fm <- new(Class='FBAmodel', S=s.al, rxnData=new("AnnotatedDataFrame", data=rxns.al, varMetadata=data.frame(labelDescription=colnames(rxns.al), row.names=colnames(rxns.al))), 
#gpr=new("AnnotatedDataFrame", data=gpr, varMetadata=data.frame(labelDescription=colnames(gpr), row.names=colnames(gpr))))

#ss.fbam(fm, r='FMN')
