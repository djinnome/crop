##jmd
##2.7.12
##which_col.r

##get column indices of 1's ordered by their row
#m is logical matrix
which.col <- function(m){
 w <- which(m, arr.ind=TRUE)
 colnames(w) <- c('row', 'col') #lose colnames of w when using Matrix pkg classes
 w.o <- w[order(w[,'row']),]
 return(w.o[,'col'])
}
