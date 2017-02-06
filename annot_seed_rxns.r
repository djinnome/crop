##jmd
##5.31.11
##annot_seed_rxns.r
##assign probs to rxns using eficaz, thermo, & nc.pwy

require('Matrix')
options(stringsAsFactors=FALSE)

setwd('/msc/neurospora/FBA/seed')

source('../farm/imbalance_seed.r')
source('../farm/probs4rxns.r')

##read
bm <- read.delim('biomass_seed.txt')
vogel <- read.delim('vogel_seed.txt')
rxns <- read.csv('/msc/neurospora/FBA/seed/ModelSEED-reactions.csv', as.is=TRUE)
rownames(rxns) <- rxns$DATABASE
rxns.meta <- read.delim('/msc/neurospora/FBA/farm_data/rxn_annot.txt', as.is=TRUE)
nc.pwy <- read.delim('/msc/neurospora/FBA/Neurospora/Nogenerics/nc10cyc.pwy')
s <- read.delim('s_seed.txt', as.is=FALSE, na='NIL')
s.sp <- sparseMatrix(i=as.numeric(s$compound), j=as.numeric(s$rxn), x=as.real(s$coeff))
dimnames(s.sp) <- list(levels(s$compound), levels(s$rxn))
nmets <- nrow(s.sp); nrxns <- ncol(s.sp)

##add meta column
#take 1st metacyc rxn that matches
rxns$meta <- NA
rxns$meta[!is.na(rxns$KEGG)] <- apply(X=rxns[!is.na(rxns$KEGG),], 1, FUN=function(x){ rxns.meta$FRAME.ID[grep(x=rxns.meta$Kegg, pattern=x['KEGG.ID.S.'])][1] })

##get eficaz probs
#ef.prob.v <- ef.probs('/msc/neurospora/FBA/Neurospora/NC10_CALLGENES_FINAL_2.ecpred', low.conf.prob=0.3)
#write.csv(round(ef.prob.v, 2), 'ec_probs.csv')
prob.mat <- read.csv('../farm_data/ec_probs.csv')
dimnames(prob.mat) <- list(prob.mat[,1], c('EC', 'prob'))
prob.mat.4l <- prob.mat[grep('.+\\..+\\..+\\..+', prob.mat$EC),]

##match EC numbers for prob & nc.pwy
rxns[rxns$EC=="", 'EC.NUMBER.S.'] <- NA
nc.pwy$ec <- rxns.meta$EC[match(nc.pwy$rxn, rxns.meta$rxn)]
rxns$prob <- rxns$nc.pwy <- 0
for (i in 1:nrow(rxns)){
    if (!is.na(rxns$EC[i])){
        #ec.tmp has some '', but that's ok
        ec.tmp <- unlist(strsplit(x=rxns$EC[i], split=',|\\|'))
        #if ec matches nc.pwy, give nc.pwy=0.5
        if (any(ec.tmp %in% nc.pwy$ec)){ rxns$nc.pwy[i] <- 0.5 }
        #get ec.prob
        pm.match <- which(apply(prob.mat.4l, MARGIN=1, FUN=function(p){ p['EC'] %in% ec.tmp }))
        rxns$prob[i] <- 1-prod(1-prob.mat.4l$prob[pm.match])
    }
}
#if meta rxn name assoc w. meta frame id is in nc.pwy, nc.pwy=1
rxns$nc.pwy[rxns.meta$rxn[match(rxns$meta, rxns.meta$FRAME, incomparables=NA)] %in% nc.pwy$rxn] <- 1
#prob(nc.pwy)=1
rxns$prob[rxns$nc.pwy==1] <- 1

##thermo probs

##match rxns to s.sp columns
rxns.match <- rxns[gsub('-L2R|-R2L', '', colnames(s.sp)),]
rownames(rxns.match) <- colnames(s.sp)
rxns.match$DATABASE[is.na(rxns.match$DATABASE)] <- rownames(rxns.match)[is.na(rxns.match$DATABASE)]
rxns.match$prob[is.na(rxns.match$prob)] <- 0
#reversible
bare.rxns <- gsub('-L2R$|-R2L$', '', rownames(rxns.match))
rxns.match$rev <- as.numeric(bare.rxns %in% bare.rxns[duplicated(bare.rxns)])
#ASSUME: gibbs is opposite sign if R2L, i.e. deltaG given for L2R independently of listed thermo feasibility
rxns.match1 <- rxns.match
rxns.match1$DELTAG[grep('-R2L$', rownames(rxns.match))] <- -1*rxns.match$DELTAG[grep('-R2L$', rownames(rxns.match))]

##unbalanced probs
#UDP -> UDPglucose had prob=1
cpds <- read.csv('ModelSEED-compounds.csv')
cpds.match <- cpds[match(rownames(s.sp), cpds$PRIMARY.NAME),]
rownames(cpds.match) <- rownames(s.sp)
#get unbalanced rxns
#met.mat <- get.unb.rxns(cpds=cpds.match)
#write.csv(met.mat, 'met_mat_seed.csv')
met.mat <- read.csv('met_mat_seed.csv')
#imbalance matrix
imb.mat <- t(s.sp) %*% met.mat
#penalize probs of unbalanced
rxns.match1$prob[rowSums(abs(imb.mat))>0] <- -1000

##special rxns
rxns.match1[c(vogel$rxn, unique(bm$rxn)), 'prob'] <- 1
  
##write
write.csv(rxns.match1, 'rxns_seed_annot.csv')

###check#########################################################################################
all(colnames(s.sp)==rownames(rxns.match1)) #TRUE
