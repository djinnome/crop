##jmd
##4.17.12
##get_gene_supp_mat.r

split.names <- function(x){ unlist(strsplit(x=gsub('[CCO-EXTRACELLULAR]', '', x, fixed=TRUE), split='\\.5\\.' )) }

multi.met.gene.supp <- function(z0){
    z <- t(z0)
    multi.met <- grep(', ', rownames(z), value=TRUE, fixed=TRUE)
    for (mm in multi.met){
        mets <- unlist(strsplit(x=mm, split=', ', fixed=TRUE))
        mets <- intersect(mets, rownames(z))
        purple.ind <- which(z[mm,]==3)
        if (length(mets)>1){
            z[mm, colSums(z[mets,]==1)>0] <- 4
            z[mm, colSums(z[mets,]==2)>0] <- 5
            #if one met predicted and other verified, then combo is predicted & verified
            z[mm, colSums(z[mets,]==1)>0 & colSums(z[mets,]==2)>0] <- 6
            z[mm, colSums(z[mets,]==3)>0] <- 6
        } else {
            z[mm, which(z[mets,]==1)] <- 4
            z[mm, which(z[mets,]==2)] <- 5
            z[mm, which(z[mets,]==3)] <- 6
        }#end else
        z[mm, purple.ind] <- 3
    }#end for
    return(t(z))
}

#3=pred&obs;2=pred only;1=obs only;0=neither
get.gene.supp.mat <- function(test.supp.df, sp, rxns.lst, supp.col='SUPPLEMENTS', obs.supp.names, rm.cpds='CPD2T-61', ub=named.vec(10**3, names=colnames(sp)), supp.ub=10, eps=10**-6){
    ##get pred
    u.supps <- gsub('\\(|\\)', '', unique(unlist(strsplit(split=' or ', x=test.supp.df$SUPP))))
    u.supps <- setdiff(u.supps, rm.cpds)    
    test.supp.df2 <- test.supp.df
    #loop thru genes & make SUPPLEMENTS = setdiff(u.supps, vector of old supps split by ' or '), all pasted together by ' or '
    test.supp.df2[,supp.col] <- apply(test.supp.df, 1, FUN=function(x){ paste(u.supps, collapse=' or ') })
    tes.nonsupp <- test.supp(sp=sp, ko.lst=rxns.lst[rownames(test.supp.df2)], annot=test.supp.df2, sense='E', ub=ub, supp.ub=supp.ub)
    pred.supps <- unlist(tes.nonsupp[[2]])[unlist(tes.nonsupp[[2]])>eps]
    
    ##get matrix of obs & prediction pairs
    pred.supp.mat <- apply(as.matrix(names(pred.supps)), 1, FUN=function(x) split.names(x))  
    pred.supp.mat[1,] <- paste(pred.supp.mat[1,], '.5', sep='')
    
    obs.supp.mat <- apply(as.matrix(obs.supp.names), 1, FUN=function(x) split.names(x))  
    obs.supp.mat[1,] <- paste(obs.supp.mat[1,], '.5', sep='')
    obs.supp.mat <- obs.supp.mat[,!(obs.supp.mat[2,] %in% rm.cpds)]
    
    ##make mat
    #3=pred&obs;2=pred only;1=obs only;0=neither
    gene.supp.mat <- matrix(0, nrow=nrow(test.supp.df), ncol=length(u.supps), dimnames=list(rownames(test.supp.df), gsub(' and ', ', ', u.supps)))
    #pred
    for (i in 1:ncol(pred.supp.mat)){ 
        gene.supp.mat[pred.supp.mat[1,i], pred.supp.mat[2,i]]  <- gene.supp.mat[pred.supp.mat[1,i], pred.supp.mat[2,i]]+2
    }
    #obs
    for (i in 1:ncol(obs.supp.mat)){ 
        gene.supp.mat[obs.supp.mat[1,i], obs.supp.mat[2,i]]  <- gene.supp.mat[obs.supp.mat[1,i], obs.supp.mat[2,i]]+1
    }
    
    gsm <- multi.met.gene.supp(gene.supp.mat)
    return(gsm)
}
