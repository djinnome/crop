##jmd
##5.31.11
##imbalance_seed.r

require('CHNOSZ')

get.unb.rxns <- function(cpds, els=c('C', 'N', 'O', 'P', 'S')){
    #metabolite element matrix
    met.mat <- matrix(0, nrow=nrow(cpds), ncol=length(els), dimnames=list(cpds$PRIMARY.NAME, els))
    for (i in 1:nrow(cpds)){
        form <- cpds[i, 'FORMULA']
        if (form!='' & !is.na(form)){
            if (length(grep('*', form, fixed=TRUE))>0){
                ss <- unlist(strsplit(form, split='*', fixed=TRUE))
                form <- paste('(', ss[1], ')', ss[2], sep='')
            }
            mm <- makeup(makeup(form),"")
            mm.match <- mm[names(mm) %in% els]
            met.mat[i, names(mm.match)] <- as.numeric(mm.match)
        }
        if (i %% 1000 == 0) print(i)
    }
    return(met.mat)
}
