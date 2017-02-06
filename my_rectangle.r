##jmd
##4.19.12
##my_rectangle.r

my.rectangle <- function(z, equal.val, border='gray', col='purple', clear.prev.sq=TRUE, density=25, ...){
    arr.inds <- which(z==equal.val, arr.ind=TRUE)
    if (nrow(arr.inds)>0){
        x <- arr.inds[,1]; y <- arr.inds[,2]
        #clear square of previous color
        if (clear.prev.sq){
            rect(x-0.5, y-0.5, x+0.5, y+0.5, col='white', border=border)
        }
        #draw in
        rect(x-0.5, y-0.5, x+0.5, y+0.5, col=col, density=density, border=border, ...)
    }#end if
}#end fcn
