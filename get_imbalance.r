##jmd
##2.7.12
##get_imbalance.r

get.imbalance <- function(sp, met.mat, biomass.rxns, trans.regexp='-TRANS-RXN'){
    bal.mat <- as.matrix(t(t(met.mat) %*% sp))
    trans.rxns <- grep(trans.regexp, colnames(sp), value=TRUE)
    bal.mat[intersect(rownames(bal.mat), c(trans.rxns, biomass.rxns)),] <- 0
    bal.mat <- bal.mat[,colSums(abs(bal.mat))>0]
}
