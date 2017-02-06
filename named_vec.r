##jmd
##2.21.12
##named_vec.r

named.vec <- function(data=NA, names){
    vec <- rep(NA, length(names))
    vec[1:length(vec)] <- data
    names(vec) <- names
    return(vec)
}
