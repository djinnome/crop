##jmd
##4.23.12
##sl_heat.r

setwd("C:/Documents and Settings/jdreyf/Desktop/farm_data")
source('C:/Documents and Settings/jdreyf/Desktop/farm/make_grid.r')
source('C:/Documents and Settings/jdreyf/Desktop/farm/my_rectangle.r')

sl <- read.csv('pred_synth_lethals.csv', as.is=TRUE)
#sl <- t(apply(sl0, 1, FUN=function(x){ x <- as.character(x); c(unlist(strsplit(x[1], split=',')), unlist(strsplit(x[2], split=','))) }))
#add pyr1;uc5
sl <- rbind(sl, c('pyr-1', 'uc-5', 'NCU06532.5', 'NCU07334.5'))
rownames(sl) <- paste(sl[,1], ',', sl[,2], sep='')

map <- data.frame(sym=c(sl[,1], sl[,2]), ncu=c(sl[,3], sl[,4]))
map <- map[!duplicated(map),]

#pwys.sl <- read.delim('pwys-of-synthetic_lethals.tsv', as.is=TRUE)
pwys.sl <- read.delim('pred_synth_lethals_pwys.tsv', as.is=TRUE)
#add pyr1;uc5
pwys.sl <- rbind(pwys.sl, c('NCU06532.5', 'NCU07334.5', 'No pathway'))
rownames(pwys.sl) <- paste(pwys.sl[,1], ',', pwys.sl[,2], sep='')
#check that pwys.sl matches sl: want these to be 0
length(setdiff(paste(sl[,3], ',', sl[,4], sep=''), rownames(pwys.sl))) #0!
length(setdiff(rownames(pwys.sl), paste(sl[,3], ',', sl[,4], sep=''))) #0!

#sl.names <- unique(unlist(strsplit(x=rownames(sl), split=',')))
sl.names <- unique(c(sl[,1], sl[,2]))
sl.mat <- matrix(0, nrow=length(sl.names), ncol=length(sl.names), dimnames=list(sl.names, sl.names))

#1: common; 2: diff; 3: not in pwys
for (i in 1:nrow(sl.mat)){
    for (j in 1:ncol(sl.mat)){
        if(paste(sl.names[i], ',', sl.names[j], sep='') %in% rownames(sl) | paste(sl.names[j], ',', sl.names[i], sep='') %in% rownames(sl)){
            ncu.pair <- map[map$sym %in% c(sl.names[i], sl.names[j]), 2]
            pwy.sl.row <- which(rownames(pwys.sl) %in% paste(ncu.pair, collapse=',') | rownames(pwys.sl) %in% paste(ncu.pair[2:1], collapse=','))
            if (pwys.sl[pwy.sl.row, 'Category']=='Common pathway') sl.mat[i,j] <- 1
            if (pwys.sl[pwy.sl.row, 'Category']=='Different pathway') sl.mat[i,j] <- 2
            if (pwys.sl[pwy.sl.row, 'Category']=='No pathway') sl.mat[i,j] <- 3
        }#end if paste
    }#end if j
}#end if i
diag(sl.mat) <- -1

hc <- hclust(dist(t(sl.mat)))
z <- sl.mat[hc$order, hc$order]
make.grid(z, color.v=c('gray', 'white', 'magenta', 'cyan', 'orange'), cex.axis=0.3, lwd=1)

z2 <- matrix(0, nrow=nrow(z), ncol=ncol(z), dimnames=dimnames(z))
z2[rbind(c('am', "en(am)-2"), c("en(am)-2", 'am'), c('pro-3', 'ota'), c('ota', 'pro-3'), c('pyr-1', 'uc-5'), c('uc-5', 'pyr-1'))] <- 1
my.rectangle(z2, equal.val=1, border='black', col=NULL, clear.prev.sq=FALSE, density=0)
