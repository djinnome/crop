##jmd
##12.5.12
##sl_add_names.r

sl.add.names <- function(sl, map, symbols.col='Symbols'){
    #get unique names from 1st two columns of sl
    sl.names <- sort(unique(c(sl[,1], sl[,2])))
    names(sl.names) <- sl.names
    names.in.map <- intersect(rownames(map), names(sl.names))
    sl.names[names.in.map] <- ncu.gene.map[names.in.map, symbols.col]
    sl.names.mat <- cbind(sl.names[sl[,1]], sl.names[sl[,2]])
    dimnames(sl.names.mat) <- list(1:nrow(sl.names.mat), paste('GeneName', 1:2, sep=''))
    return(sl.names.mat)
}
