##jmd
##10.7.11
##parse_supp.r

parse.supp <- function(annot){
    lst2 <- list()
    for (j in 1:nrow(annot)){
        if (annot$SUPPLEMENT[j]!=''){
            y <- strsplit(annot$SUPPLEMENT[j], split=' or ')[[1]]
            z <- list()
            for (i in 1:length(y)){
                z[[i]] <- unlist(strsplit(gsub('\\(|\\)', '', y[i]), split=' and '))
            }
            lst2[[ annot[j, 'Locus'] ]] <- z
        }#end if
    }#end for
    return(lst2)
}
