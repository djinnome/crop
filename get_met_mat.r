##jmd
##1.30.12
##get_met_mat.r

#smm.match matches smm to s.rownames
#this fcn is for metacyc only right now. in future, this can call get.mm.meta or get.mm.seed, depending on 'db' parameter
#colsums.mm.thresh represents the min number of mets an element needs to be involved in
get.met.mat <- function(univ.smm, org.smm, s.rownames, chem.form.col='CHEMICAL.FORMULA', colsums.mm.thresh=2, suffix.regexp='\\[CCO-.+\\]'){    
    ##combine
    smm <- merge.smm(univ.smm, org.smm)
    
    ##match smm for normal & compartmentalized mets, too
    smm.match <- smm[gsub(suffix.regexp, '', s.rownames),]
    rownames(smm.match) <- smm.match$FRAME <- s.rownames
    
    ##parse
    smm.match <- clean.meta.chem.form(smm.match)
    #check
    #cat('proportion of metabolites w/o chemical formula:', mean(is.na(smm.match$CHEM)), '\n')
    #get unique elements
    els.w.count <- chemForm2df(smm.match[,chem.form.col])
    u.els <- els.w.count$els
    
    ##make matrix with metabs as rows and molecules as columns
    #can sum rows according to s to check balance
    #using matrix() instead of Matrix makes it much faster
    met.mat <- matrix(0, nrow=nrow(s.sp), ncol=length(u.els), dimnames=list(s.rownames, u.els))
    for (i in 1:nrow(met.mat)){
        if (smm.match[i, chem.form.col]!=''){
            df.tmp <- chemForm2df(smm.match[i, chem.form.col])
            met.mat[i, df.tmp$els] <- df.tmp$count
        }
    }
    
    ##subset & order
    mm.ss <- met.mat[,colSums(met.mat)>=colsums.mm.thresh]
    mm.ss <- as.matrix(mm.ss[,order(-colSums(mm.ss))])
    
    return(mm.ss)
}
