##jmd
##3.3.11
##probs4rxns.r

ef.probs <- function(ecpred.path, low.conf.prob=0.3){
    ef0 <- read.delim(ecpred.path, colClasses=c(rep('character', 6), 'numeric', 'numeric'))
    ef <- ef0[ef0$ec_number!='',]
    ef <- ef[order(ef$ec_num),]
    #LOW conf gets 30%
    ef$mean[ef$conf=='LOW'] <- low.conf.prob
    #combine probs
    probs <- tapply(ef$mean, INDEX=ef$ec_num, FUN=function(x){ 1-prod(1-x) })
    prob.mat <- data.frame(EC=rownames(probs), prob=probs)
    dimnames(prob.mat) <- list(prob.mat[,1], c('EC', 'ef.prob'))
    return(prob.mat)
}

#currently HARD-CODES ef.probs & thermo.probs
#still need to decide how to structure this feature vector
#don't use pwy yet & gpr seems redundant w/ eficaz
#gpr.prob = P(rxn in gpr|rxn in org); pwy.prob = P(rxn in pwy|rxn in org)
get.probs <- function(rxn.annot, prob.cols=c('ef.prob', 'thermo.prob'), gpr.prob=0.8, pwy.prob=0.9, unb.rm.pen=1, h.unb.add.pen=10, unb.add.pen=10**3){
    gpr.coeff <- gpr.prob*rxn.annot$gpr+(1-gpr.prob)*(1-rxn.annot$gpr)
    
    rxn.annot$ef.prob[which(is.na(rxn.annot$ef.prob))] <- 0
    rxn.annot$thermo.prob[which(is.na(rxn.annot$thermo.prob))] <- mean(rxn.annot$thermo.prob, na.rm=TRUE)
    
    probs <- apply(rxn.annot[,prob.cols], MARGIN=1, FUN=prod)
    #probs <- probs*gpr.coeff
    #special cases
    probs[rxn.annot$unb.rm>0 | rxn.annot$h.unb<0] <- probs[rxn.annot$unb.rm>0 | rxn.annot$h.unb<0]-unb.rm.pen
    probs[rxn.annot$h.unb>0] <- probs[rxn.annot$h.unb>0]-h.unb.add.pen
    probs[rxn.annot$unb.add>0] <- -unb.add.pen
    
    #probs[is.na(probs)] <- 0
    
    return(probs)
}



    ##get probs
    #special rxns
    #rxns$prob[rxns$rxn %in% c(bm$rxn, bm.goals$rxn, met.rxns, vogel$rxn)] <- 1
    
    #na probs
    #rxns$prob[is.na(rxns$prob)] <- 0
    #rxns in gpr but w/o eficaz
    #rxns$prob[rxns$rxn %in% gpr$RXN & rxns$prob==0] <- 0.25
    #prob(nc.pwy)=x+(1-x)*max(prob, 0), where max is for unbalanced nc.pwy rxns w/ prob<0 
    #rxns$prob[rxns$nc.pwy==1 & rxns$prob>=0] <- 0.5+0.5*apply(as.matrix(rxns$prob[rxns$nc.pwy==1 & rxns$prob>=0]), 1, FUN=function(x){ max(x, 0) })
