##jmd
##4.18.12
##make_grid.r

make.grid <- function(z, mar=c(9,5,4,2), color.v=c('white', 'red', 'blue', 'purple'), cex.axis=0.6, lwd=0, lwd.ticks=0){
    par(mar=mar+0.1)
    image(x=1:nrow(z), y=1:ncol(z), z=z, col=color.v, xlab=NA, ylab=NA, axes=FALSE)
    abline(h=2:ncol(z)-0.5, col=c('darkgray', 'gray'))
    abline(v=2:nrow(z)-0.5, col=c('darkgray', 'gray'))
    axis(1, at=0:(nrow(z)+1), labels=c('', rownames(z), ''), las=2, cex.axis=cex.axis, tck=0, lwd=lwd, lwd.ticks=lwd.ticks)
    axis(2, at=0:(ncol(z)+1), labels=c('', colnames(z), ''), las=2, cex.axis=cex.axis, tck=0, lwd=lwd, lwd.ticks=lwd.ticks)
    axis(3, at=0:(nrow(z)+1), labels=rep('', nrow(z)+2), tck=0, lwd=lwd, lwd.ticks=lwd.ticks)
    axis(4, at=0:(ncol(z)+1), labels=rep('', ncol(z)+2), tck=0, lwd=lwd, lwd.ticks=lwd.ticks)
}
