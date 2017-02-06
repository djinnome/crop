##jmd
##3.5.12
##biolog_confusion.r

#na.rm asks whether or not to remove NAs for computation. if TRUE, then growth in 1 experiment & NA in another gives overall 'growth' else it gives overall NA.
biolog.confusion <- function(growth.low=2, nogrowth.high='<1', na.rm=TRUE,
exp.biolog.file=c('pm1-vogels-0.1-day4.tsv', 'pm1-water-0.01-day3.tsv'), pred.biolog.file='biolog_predict.tsv', 
levels.ordered=c(0, '<1', 1, '1+', 2, '2+', 3, '3+')){
    
    dat.lst <- list()
    bin.growth.mat <- NULL
    for (i in 1:length(exp.biolog.file)){
        dat.lst[[i]] <- read.delim(exp.biolog.file[i])
        colnames(dat.lst[[i]]) <- c('well', 'medium', 'growth')
        rownames(dat.lst[[i]]) <- dat.lst[[i]]$well
        dat.lst[[i]]$growth <- factor(dat.lst[[i]]$growth, levels=levels.ordered, ordered=TRUE)
        
        dat.lst[[i]]$bin.growth <- NA
        dat.lst[[i]]$bin.growth[dat.lst[[i]]$growth<=nogrowth.high] <- 0
        dat.lst[[i]]$bin.growth[dat.lst[[i]]$growth>=growth.low] <- 1
        
        bin.growth.mat <- cbind(bin.growth.mat, dat.lst[[i]]$bin.growth)
    }
    dat <- dat.lst[[1]]
    dat$growth <- apply(bin.growth.mat, 1, FUN=function(x){
        condn <- all(x==mean(x, na.rm=na.rm), na.rm=na.rm)
        if (!is.na(condn) & condn){ mean(x, na.rm=TRUE) } else { NA }
    })
        
#    g <- named.vec(NA, rownames(dat.lst[[1]]))
#    g[dat$growth<=nogrowth.high] <- 0
#    g[dat$growth>=growth.low] <- 1
    
    pred <- read.delim(pred.biolog.file)
    #make sure only have wells for the same plate as experiment
    rownames(pred) <- pred$well

    int <- intersect(pred$well, rownames(dat))
    mat <- data.frame(cpd=substr(x=pred[int, 'cpd'], start=1, stop=50), pred=pred[int, 'pred'], obs=dat[int, 'growth'])
    tab <- table(pred=pred[int, 'pred'], obs=dat[int, 'growth'])
    tab.out <- rbind(tab, c(tab[1,1]/sum(tab[,1]), tab[2,2]/sum(tab[,2])))
    return(list(tab=tab.out, mat=mat))
}
