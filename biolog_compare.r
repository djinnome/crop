##jmd
##4.26.12
##biolog_compare.r

source('/msc/neurospora/FBA/farm/farm_header.r')

EPS <- 0.02

##get pred
##run run.pred.biolog
bb <- read.delim('biolog_non-imported_model_cpds.tsv', as.is=TRUE)

pb <- run.pred.biolog(s.al, default.nuts=vogel, add.trans.mets=bb$cpd[bb$should.be.imported=='yes'], 
na.is.nogrowth=TRUE, extracellular=TRUE, write.out.tsv='biolog_predict_fba', eps=EPS); pb[1:4,]

#pred <- read.delim('biolog_predict_fba.tsv')
pred <- pb
#check quant.pred
#sort(pred$quant.pred[!is.na(pred$quant.pred) & pred$quant.pred>0]) #none near 10**-6, max is 9
pred$pred[is.na(pred$quant.pred)|pred$quant.pred<EPS] <- 0
rownames(pred) <- paste(pred$plate, pred$well, sep='.')

##get matrix of obs
exp.mat.all <- NULL
for (plate in paste('PM', 1:4, sep='')){
    ##get matrix of scores for matching PM plate
    exp.tmp <- read.csv(paste(tolower(plate), '_24april12.csv', sep=''))
    rownames(exp.tmp) <- paste(plate, gsub('(0)(.)', replacement='\\2', exp.tmp$well), sep='.')
    exp.mat <- as.matrix(exp.tmp[,grep('rep', colnames(exp.tmp), ignore.case=TRUE)])
    #translate '<' and '+'
    for (i in 1:ncol(exp.mat)){
        minus.ind <- grep('<', exp.mat[,i], fixed=TRUE)
        plus.ind <- grep('+', exp.mat[,i], fixed=TRUE)
        
        exp.mat[minus.ind, i] <- as.numeric(substr(x=exp.mat[minus.ind, i], start=2, stop=2))-0.25
        exp.mat[plus.ind, i] <- as.numeric(substr(x=exp.mat[plus.ind, i], start=1, stop=1))+0.25
    }
    exp.mat <- matrix(as.numeric(exp.mat), ncol=4, dimnames=list(rownames(exp.mat), paste('rep', 1:4, sep='')))
    exp.mat.all <- rbind(exp.mat.all, exp.mat)
}
#add '0' to 1-digit well names
rownames(exp.mat.all) <- gsub('(PM...)(.)$', replacement='\\10\\2', rownames(exp.mat.all))
 
##summarize obs scores per row
growth.low=2; nogrowth.high <- 0.9; nogrowth.low <- 0
exp.growth <- rowMeans(exp.mat.all)
exp.growth.bin <- rep(NA, length(exp.growth))
names(exp.growth.bin) <- names(exp.growth)
exp.growth.bin[exp.growth>=growth.low] <- 1
exp.growth.bin[nogrowth.low<=exp.growth & exp.growth<=nogrowth.high] <- 0
summary(as.factor(exp.growth.bin[names(exp.growth.bin) %in% rownames(pred)]))

##compare
int <- intersect(names(exp.growth), rownames(pred))
mat <- data.frame(cpd=substr(x=pred[int, 'cpd'], start=1, stop=50), pred=pred[int, 'pred'], obs=exp.growth.bin[int])
tab <- table(pred=mat[,'pred'], obs=mat[,'obs'])
tab.out <- rbind(tab, c(tab[1,1]/sum(tab[,1]), tab[2,2]/sum(tab[,2])))
print(tab.out)

##check
mm <- mat[rowSums(is.na(mat))==0,]
mm <- mm[mm$pred!=mm$obs,]
wrong.df <- data.frame(cpd=mm$cpd, pred=pred[rownames(mm), 'pred'], quant.pred=signif(pred[rownames(mm), 'quant.pred'], 2), obs=mm$obs, exp.mat.all[rownames(mm),])
wrong.df <- wrong.df[order(wrong.df$pred, rownames(wrong.df)),]
wrong.df
