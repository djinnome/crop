##jmd
##5.25.11
##subset_mat.r
##fcn to 'see' non-zero subset of sparse matrix

ss.mat <- function(sp, c=NULL, r=NULL){
    if(is.null(r)){ r.ss <- 1:nrow(sp) } else { r.ss <- r }
    if(is.null(c)){ c.ss <- 1:ncol(sp) } else { c.ss <- c }
    ss <- sp[r.ss, c.ss]
    if (length(c)==1|length(r)==1){ ss <- as.matrix(ss) }
    ss <- ss[rowSums(abs(ss))>0, colSums(abs(ss))>0]
    return(ss)
}

##get matrix SubSet of Rxns Of a Metabolite
ss.rom <- function(sp, met.v){ 
    ss.mat(sp=sp, c=rownames(as.matrix(ss.mat(sp=sp, r=met.v)))) 
}
