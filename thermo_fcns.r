##jmd
##8.4.11
##thermo_fcns.r

#need to strip compartmentalization from rownames(s)
get.dg0 <- function(dg0.path="/msc/neurospora/FBA/Meta/Gibbs.lisp", sp){
  ##parse thermo metabolites (jeremy left off direction in metacyc-deltag0.txt)
  dg0.mets <- as.matrix(read.delim(dg0.path))
  dg0.met.mat <- t(apply(dg0.mets, 1, FUN=function(x){ 
    unlist(strsplit(x=gsub('\\(modify-cpd \'|:estimated-energy |:estimated-uncertainty | :citations.+','', x), split=" ")) 
  }))
  #fix quotes here: '
  #get rid of pipes in met names, eg |Pi|
  dg0.met.mat[,1] <- gsub('\\|', '', dg0.met.mat[,1])
  dimnames(dg0.met.mat) <- list(dg0.met.mat[,1], c('mol', 'dg0', 'sd.dg0'))
  dmm <- as.numeric(dg0.met.mat[,'dg0'])
  names(dmm) <- rownames(dg0.met.mat)
  dmm <- dmm[gsub('\\[CCO-.+\\]$', '', rownames(sp))]
  #treat dg0=-9.53 as NA for molecules other than 'PROTON', since this value appears for molecules w/ empty mol file
  dmm[dmm==-9.53 & names(dmm)!='PROTON'] <- NA
  #assume mets w/o gibbs show up on both sides and approximately cancel, so treat their gibbs as 0
  #dmm[is.na(dmm)] <- 0  
  
  ##compute dg0 for rxns
  #only have dg for 2/3 of metacyc mets, so will have many rxns w/o dg
  dg0.rxns <- as.numeric(t(sp) %*% dmm)
  names(dg0.rxns) <- colnames(sp)
  return(dg0.rxns)
}

##obtain rxn (directional) probabilities from dg0
#use sd=10.25 since 20.5kJ/mol corresponds to Jankowski estimate of 2.2kcal/mol, which represents 1 std err ~= 2 sd 
#although i think c henry says in his email that he would prefer that this be 1 sd
#problem: using uncertainty removes almost all thermo penalties, so use sd=10.25/4 here (arbitrarily)
#known.mets.conc is vector w/ names %in% rownames(sp) and entries = molarity; w/ proton, it changes E(thermo.prob) by ~0.1 & adds ~1000 rxns to those w/ p<0.1
dg2prob <- function(sp, dg0.rxns, known.mets.conc=NULL, dg0.sd=10.25/4, use.met.conc.sd=FALSE){
  #if know some concentrations
  if (!is.null(known.mets.conc)){ 
    dg.known.mean <- 2.5*as.numeric(t(sp[names(known.mets.conc),]) %*% log(known.mets.conc)) 
    unknown.mets <- setdiff(rownames(sp), names(known.mets.conc))
  } else { 
    dg.known.mean <- numeric(ncol(sp))
    names(dg.known.mean) <- colnames(sp)
    unknown.mets <- rownames(sp)
  }#end else

  rs.neg <- apply(sp[unknown.mets,], 2, FUN=function(x){ sum(abs(x[x<0])) })
  rs.pos <- apply(sp[unknown.mets,], 2, FUN=function(x){ sum(x[x>0]) })
  #vector operation
  dg.mean <- dg.known.mean+dg0.rxns+2.5*(3.5*rs.neg-13.8*rs.pos)
  if (!use.met.conc.sd){ dg.sd <- dg0.sd } else { dg.sd <- sqrt(dg0.sd^2+0.04*((rs.pos-rs.neg)^2+2*(rs.pos+rs.neg)^2)) }
  #thermo prob
  thermo.prob <- pnorm(0, mean=dg.mean, sd=dg.sd)

  return(thermo.prob)
}
