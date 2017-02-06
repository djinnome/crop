##jmd
##5.27.11
##sparse_mat_fcns.r

##regular matrix to sparse format
dense2sp <- function(m, rxn.name='rxn', cpd.name='cpd', coeff.name='coeff'){
    w <- which(m!=0, arr.ind=TRUE)
    sp <- data.frame(colnames(m)[w[,2]], rownames(m)[w[,1]], m[w])
    colnames(sp) <- c(rxn.name, cpd.name, coeff.name)
    return(sp)
}

##half-eqn to sparse format
halfEqn2sp <- function(a, sgn=1){
    els <- gsub('^\\([0123456789.]+\\) ', '', a)
    coeff <- rep(1, length(a))
    non.one <- grep('^\\([0123456789.]+\\) ', a)
    if (length(non.one)>0){ coeff[non.one] <- apply(X=as.matrix(a[non.one]), MARGIN=1, FUN=function(x){ gsub('\\(|\\)', '', unlist(strsplit(x, split=' '))[1]) }) }
    s.sp.tmp <- data.frame(compound=els, coeff=sgn*as.numeric(coeff))
    if (!all(is.na(a))) return(s.sp.tmp) else return(NULL)
}

##turn chris henry's seed spreadsheet into sparse matrix format
#did some manual work: filled in empty THERM to '<=>' to be safe, & turned some '<=' or '=>' into '<=>' in 'NAME.EQ'
spread2sparse <- function(spr){
    s.sp <- NULL
    for (i in 1:nrow(spr)){
        #have to use 'eqn' instead of 'name.eq' since some name.eq='NONE'
        ss <- unlist(strsplit(x=spr$EQUATION[i], split=' <=> ', fixed=TRUE))
        a <- unlist(strsplit(x=ss[1], split=' + ', fixed=TRUE))
        b <- unlist(strsplit(x=ss[2], split=' + ', fixed=TRUE))
        if (spr$THERM[i] %in% c('<>', '>')) s.sp <- rbind(s.sp, cbind(rxn=paste(spr$DATABASE[i], '-L2R',sep=''), rbind(halfEqn2sp(a, sgn=-1), halfEqn2sp(b))))
        if (spr$THERM[i] %in% c('<>', '<')) s.sp <- rbind(s.sp, cbind(rxn=paste(spr$DATABASE[i], '-R2L',sep=''), rbind(halfEqn2sp(a), halfEqn2sp(b, sgn=-1))))
        if (i %% 1000 == 0) print(i)
    }
    return(s.sp)
}
